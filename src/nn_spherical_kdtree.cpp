// file: nn_spherical_kdtree.cpp

#include <Rcpp.h>
#include <nanoflann.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace nanoflann;


// ============================================================================
// Point cloud
// ============================================================================

struct PointCloud3D {

  const NumericVector& x;
  const NumericVector& y;
  const NumericVector& z;

  PointCloud3D(const NumericVector& x_,
               const NumericVector& y_,
               const NumericVector& z_)
    : x(x_), y(y_), z(z_) {}

  inline size_t kdtree_get_point_count() const {
    return x.size();
  }

  inline double kdtree_get_pt(const size_t idx, int dim) const {

    if (dim == 0) return x[idx];
    if (dim == 1) return y[idx];

    return z[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const {
    return false;
  }
};


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List nn_spherical_cpp(NumericVector sp_x,
                      NumericVector sp_y,
                      NumericVector sp_z,
                      NumericVector qp_x,
                      NumericVector qp_y,
                      NumericVector qp_z,
                      int k = 1,
                      Nullable<double> search_radius = R_NilValue,
                      int cores = 1) {
  const int n_sp = sp_x.size();
  const int n_qp = qp_x.size();
  if (sp_y.size() != n_sp || sp_z.size() != n_sp)
    stop("Sampling point vectors must have same length.");
  if (qp_y.size() != n_qp || qp_z.size() != n_qp)
    stop("Query point vectors must have same length.");

  // k == 0 is the internal representation of the R interface's k = NULL.
  const bool use_radius = search_radius.isNotNull();
  if (k < 0 || (k == 0 && !use_radius))
    stop("Supply a positive k or a search radius.");
  double radius2 = 0;
  if (use_radius) {
    const double angle = as<double>(search_radius);
    if (!std::isfinite(angle) || angle <= 0 || angle > std::acos(-1.0))
      stop("Angular radius must be in (0, pi].");
    const double chord = 2.0 * std::sin(angle / 2.0);
    radius2 = chord * chord;
    // Include the boundary despite rounding in unit-vector coordinates.
    const double eps = std::numeric_limits<double>::epsilon();
    radius2 += 16.0 * eps * radius2 + 16.0 * eps * eps;
  }
  const int limit = std::min(k, n_sp);
  std::vector<std::vector<int>> indices(n_qp);
  std::vector<std::vector<double>> distances(n_qp);

  if (n_sp > 0 && n_qp > 0) {
    PointCloud3D cloud(sp_x, sp_y, sp_z);
    typedef KDTreeSingleIndexAdaptor<
      L2_Simple_Adaptor<double, PointCloud3D>, PointCloud3D, 3, size_t> tree_t;
    tree_t tree(3, cloud, KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores > 0 ? cores : 1)
#endif
    for (int i = 0; i < n_qp; ++i) {
      const double query[3] = {qp_x[i], qp_y[i], qp_z[i]};
      std::vector<ResultItem<size_t, double>> hits;
      if (use_radius && (k == 0 || k >= n_sp)) {
        tree.radiusSearch(query, radius2, hits, SearchParameters());
      } else {
        std::vector<size_t> positions(limit);
        std::vector<double> squared(limit);
        const size_t count = tree.knnSearch(query, limit,
                                           positions.data(), squared.data());
        for (size_t j = 0; j < count; ++j) {
          if (!use_radius || squared[j] <= radius2)
            hits.emplace_back(positions[j], squared[j]);
        }
      }
      std::sort(hits.begin(), hits.end(),
        [](const ResultItem<size_t, double>& a,
           const ResultItem<size_t, double>& b) {
          return a.second == b.second ? a.first < b.first : a.second < b.second;
        });
      indices[i].reserve(hits.size());
      distances[i].reserve(hits.size());
      for (const auto& hit : hits) {
        indices[i].push_back(static_cast<int>(hit.first) + 1);
        const double half_chord = std::min(1.0, std::sqrt(hit.second) / 2.0);
        distances[i].push_back(2.0 * std::asin(half_chord));
      }
    }
  }

  // Construct R objects on the main thread, after all searches have finished.
  List output(n_qp);
  for (int i = 0; i < n_qp; ++i)
    output[i] = List::create(Named("index") = wrap(indices[i]),
                            Named("distance") = wrap(distances[i]));
  return output;
}
