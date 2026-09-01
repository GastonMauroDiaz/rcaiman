// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <cstdint>

inline uint64_t make_key(int a, int b)
{
  return (static_cast<uint64_t>(a) << 32) |
    static_cast<uint32_t>(b);
}

using namespace Rcpp;

struct Candidate {
  int from;
  int to;
  double h;
};

inline double compute_union_heterogeneity(
    double w1,
    double w2,
    double mu1,
    double mu2,
    double ssd1,
    double ssd2
) {

  return ssd1 +
    ssd2 +
    ((w1 * w2) / (w1 + w2)) *
    (mu1 - mu2) * (mu1 - mu2);
}

inline int find_root(
    std::vector<int>& parent,
    int x
) {

  while (parent[x] != x) {

    parent[x] = parent[parent[x]];

    x = parent[x];
  }

  return x;
}

inline void add_neighbor(
    std::vector<int>& v,
    int x
) {

  for (int y : v) {
    if (y == x) {
      return;
    }
  }

  v.push_back(x);
}

// [[Rcpp::export]]
IntegerMatrix merge_segments_cpp(
    NumericMatrix x_r,
    NumericMatrix x_w,
    IntegerMatrix x_seg,
    DataFrame adj,
    double scale_parameter,
    double boundary_weight,
    int max_iterations = 100000,
    bool verbose = false
) {

  const int nr = x_seg.nrow();
  const int nc = x_seg.ncol();

  // ============================================================
  // Find max id
  // ============================================================

  int max_id = 0;

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {

      int id = x_seg(i, j);

      if (id > max_id) {
        max_id = id;
      }
    }
  }

  // ============================================================
  // Region statistics
  // ============================================================

  std::vector<double> region_boundary(max_id + 1, 0.0);
  std::vector<double> region_weight(max_id + 1, 0.0);
  std::vector<double> region_mean(max_id + 1, 0.0);
  std::vector<double> region_ssd(max_id + 1, 0.0);

  std::vector<double> weighted_sum(max_id + 1, 0.0);

  std::vector<bool> active(max_id + 1, false);

  std::vector<int> parent(max_id + 1);

  for (int id = 0; id <= max_id; ++id) {
    parent[id] = id;
  }

  // ============================================================
  // First pass
  // ============================================================

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {

      int id = x_seg(i, j);

      if (id <= 0) {
        continue;
      }

      active[id] = true;

      double value = x_r(i, j);
      double weight = x_w(i, j);

      region_weight[id] += weight;
      weighted_sum[id] += weight * value;
    }
  }

  for (int id = 1; id <= max_id; ++id) {

    if (!active[id]) {
      continue;
    }

    region_mean[id] =
      weighted_sum[id] / region_weight[id];
  }

  // ============================================================
  // Second pass
  // ============================================================

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {

      int id = x_seg(i, j);

      if (id <= 0) {
        continue;
      }

      double value = x_r(i, j);
      double weight = x_w(i, j);

      double d = value - region_mean[id];

      region_ssd[id] += weight * d * d;
    }
  }

  // ============================================================
  // Adjacency graph
  // ============================================================

  IntegerVector adj_from = adj["from"];
  IntegerVector adj_to   = adj["to"];
  NumericVector adj_boundary = adj["boundary_length"];

  std::vector< std::vector<int> > neighbors(max_id + 1);

  std::unordered_map<uint64_t,double> shared_boundary;

  for (int i = 0; i < adj.nrows(); ++i) {

    int a = adj_from[i];
    int b = adj_to[i];

    add_neighbor(neighbors[a], b);
    add_neighbor(neighbors[b], a);

    uint64_t key =
      make_key(
        std::min(a,b),
        std::max(a,b)
      );

    shared_boundary[key] =
      adj_boundary[i];
  }

  // ============================================================
  // Scale normalization (Fixed)
  // ============================================================

  std::vector<double> initial_h(adj.nrows());

  for (int i = 0; i < adj.nrows(); ++i) {

    int id1 = adj_from[i];
    int id2 = adj_to[i];
    double border = adj_boundary[i];

    if (border <= 0) {
      stop("invalid border in scale normalization");
    }

    double heterogeneity = compute_union_heterogeneity(
      region_weight[id1],
                   region_weight[id2],
                                region_mean[id1],
                                           region_mean[id2],
                                                      region_ssd[id1],
                                                                region_ssd[id2]
    );

    initial_h[i] = heterogeneity / std::pow(border, boundary_weight);
  }

  std::nth_element(
    initial_h.begin(),
    initial_h.begin() + initial_h.size() / 2,
    initial_h.end()
  );

  double scale_unit =
    initial_h[initial_h.size() / 2];

  double scale_parameter_internal =
    scale_parameter * scale_unit;

  // ============================================================
  // Reusable buffers
  // ============================================================

  std::vector<int> ids;
  ids.reserve(max_id);

  std::vector<int> best_to(max_id + 1);

  std::vector<double> best_h(max_id + 1);

  std::vector<bool> used(max_id + 1);

  std::vector<Candidate> candidates;
  candidates.reserve(max_id);

  // ============================================================
  // Main loop
  // ============================================================

  int iteration = 1;

  while (iteration <= max_iterations) {

    ids.clear();

    for (int id = 1; id <= max_id; ++id) {

      if (active[id]) {
        ids.push_back(id);
      }
    }

    if (ids.size() <= 1) {
      break;
    }

    std::fill(
      best_to.begin(),
      best_to.end(),
      0
    );

    std::fill(
      best_h.begin(),
      best_h.end(),
      INFINITY
    );

    // ----------------------------------------------------------
    // Best neighbors
    // ----------------------------------------------------------

    for (int id : ids) {

      double w1 = region_weight[id];
      double mu1 = region_mean[id];
      double ssd1 = region_ssd[id];

      std::vector<int>& neighs =
        neighbors[id];

      for (int nb : neighs) {

        nb = find_root(parent, nb);

        if (nb == id)
          continue;

        if (!active[nb])
          continue;

        uint64_t key =
          make_key(
            std::min(id, nb),
            std::max(id, nb)
          );

        auto it =
          shared_boundary.find(key);

        if (it == shared_boundary.end()) {

          Rcout
          << "missing boundary "
          << id
          << " "
          << nb
          << "\n";

          stop("missing boundary");
        }

        double border = it->second;

        if (border <= 0) {
          stop("invalid border");
        }

        nb = find_root(parent, nb);

        if (nb == id) {
          continue;
        }

        if (!active[nb]) {
          continue;
        }

        double heterogeneity =
          compute_union_heterogeneity(
            w1,
            region_weight[nb],
                         mu1,
                         region_mean[nb],
                                    ssd1,
                                    region_ssd[nb]
          );

        double h = heterogeneity / std::pow(border, boundary_weight);

        if (h < best_h[id]) {

          best_h[id] = h;
          best_to[id] = nb;
        }
      }
    }

    // ----------------------------------------------------------
    // Mutual candidates
    // ----------------------------------------------------------

    candidates.clear();

    for (int id : ids) {

      int nb = best_to[id];

      if (nb == 0) {
        continue;
      }

      if (best_to[nb] != id) {
        continue;
      }

      double h = best_h[id];

      if (h > scale_parameter_internal) {
        continue;
      }

      Candidate c;

      c.from = std::min(id, nb);
      c.to = std::max(id, nb);
      c.h = h;

      candidates.push_back(c);
    }

    if (candidates.empty()) {
      break;
    }

    std::sort(
      candidates.begin(),
      candidates.end(),
      [](const Candidate& a,
         const Candidate& b) {

        return a.h < b.h;
      }
    );

    std::fill(
      used.begin(),
      used.end(),
      false
    );

    // ----------------------------------------------------------
    // Execute merges
    // ----------------------------------------------------------

    for (const Candidate& c : candidates) {

      int from = find_root(parent, c.from);
      int to = find_root(parent, c.to);

      if (from == to) {
        continue;
      }

      if (used[from] || used[to]) {
        continue;
      }

      used[from] = true;
      used[to] = true;

      double w1 = region_weight[from];
      double w2 = region_weight[to];

      double mu1 = region_mean[from];
      double mu2 = region_mean[to];

      double ssd1 = region_ssd[from];
      double ssd2 = region_ssd[to];

      double weight = w1 + w2;

      double mean_value =
        (w1 * mu1 + w2 * mu2) / weight;

      double ssd =
        compute_union_heterogeneity(
          w1,
          w2,
          mu1,
          mu2,
          ssd1,
          ssd2
        );

      region_weight[from] = weight;
      region_mean[from] = mean_value;
      region_ssd[from] = ssd;

      active[to] = false;

      parent[to] = from;

      // --------------------------------------------------------
      // Merge adjacency
      // --------------------------------------------------------

      std::vector<int>& neigh_to =
        neighbors[to];

      std::vector<int>& neigh_from =
        neighbors[from];

      uint64_t internal_key =
        make_key(
          std::min(from,to),
          std::max(from,to)
        );

      shared_boundary.erase(internal_key);

      for (int nb : neigh_to) {

        nb = find_root(parent, nb);

        if (nb == from)
          continue;

        uint64_t key_to =
          make_key(
            std::min(to, nb),
            std::max(to, nb)
          );

        double border_to = 0.0;

        auto it_to =
          shared_boundary.find(key_to);

        if (it_to != shared_boundary.end()) {

          border_to = it_to->second;

          shared_boundary.erase(it_to);
        }

        uint64_t key_from =
          make_key(
            std::min(from, nb),
            std::max(from, nb)
          );

        shared_boundary[key_from] += border_to;

        add_neighbor(neigh_from, nb);
        add_neighbor(neighbors[nb], from);
      }

      neigh_to.clear();
    }

    iteration++;
  }

  // ============================================================
  // Final relabel
  // ============================================================

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {

      int id = x_seg(i, j);

      if (id <= 0) {
        continue;
      }

      x_seg(i, j) =
        find_root(parent, id);
    }
  }

  return x_seg;
}
