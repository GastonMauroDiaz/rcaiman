// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include <unordered_map>
#include <cstdint>

using namespace Rcpp;

// ============================================================
// Edge hash
// ============================================================


inline uint64_t make_key(int a, int b)
{
  // Primero a unsigned de 32 bits, luego se expande a 64 bits de forma segura
  uint64_t ua = static_cast<uint32_t>(a);
  uint64_t ub = static_cast<uint32_t>(b);
  return (ua << 32) | ub;
}

// ============================================================
// compute_adjacency_cpp
// ============================================================

// [[Rcpp::export]]
DataFrame compute_adjacency_cpp(
    IntegerMatrix x,
    NumericMatrix w,
    bool directed = false,
    bool ignore_label_0 = true,
    bool return_border_length = true
){

  const int nr = x.nrow();
  const int nc = x.ncol();

  if (w.nrow() != nr || w.ncol() != nc)
    stop("`x` and `w` must have identical dimensions.");


  std::unordered_map<uint64_t, double> edge_lengths;

  // ==========================================================
  // Horizontal neighbors
  // ==========================================================

  for (int i = 0; i < nr; ++i) {

    for (int j = 0; j < nc - 1; ++j) {

      int a = x(i, j);
      int b = x(i, j + 1);

      // if (IntegerVector::is_na(a) ||
      //     IntegerVector::is_na(b))
      //   continue;

      if (ignore_label_0 &&
          (a == 0 || b == 0))
        continue;

      if (a == b)
        continue;

      if (!directed && a > b) {
        std::swap(a, b);
      }

      if (IntegerVector::is_na(a) ||
          IntegerVector::is_na(b))
      {
        stop("NA reached make_key()");
      }

      if (a < 0 || b < 0) {
        Rcout << "negative label: "
              << a << " "
              << b << "\n";
      }

      uint64_t key = make_key(a, b);

      double edge =
        0.5 * (
            std::sqrt(w(i,j)) +
              std::sqrt(w(i,j+1))
        );

      edge_lengths[key] += edge;

    }
  }

  // ==========================================================
  // Vertical neighbors
  // ==========================================================

  for (int i = 0; i < nr - 1; ++i) {

    for (int j = 0; j < nc; ++j) {

      int a = x(i, j);
      int b = x(i + 1, j);

      // if (IntegerVector::is_na(a) ||
      //     IntegerVector::is_na(b))
      //   continue;

      if (ignore_label_0 &&
          (a == 0 || b == 0))
        continue;

      if (a == b)
        continue;

      if (!directed && a > b) {
        std::swap(a, b);
      }

      if (IntegerVector::is_na(a) ||
          IntegerVector::is_na(b))
      {
        stop("NA reached make_key()");
      }

      if (a < 0 || b < 0) {
        Rcout << "negative label: "
              << a << " "
              << b << "\n";
      }

      uint64_t key = make_key(a, b);

      double edge =
        0.5 * (
            std::sqrt(w(i,j)) +
              std::sqrt(w(i+1,j))
        );

      edge_lengths[key] += edge;
    }
  }

  // ==========================================================
  // Output
  // ==========================================================

  const int n = edge_lengths.size();

  IntegerVector from(n);
  IntegerVector to(n);

  if (return_border_length) {

    NumericVector border(n);

    int k = 0;

    for (const auto& kv : edge_lengths) {

      uint64_t key = kv.first;

      from[k] = static_cast<int>(key >> 32);
      to[k] = static_cast<int>(key & 0xffffffff);

      border[k] = kv.second;

      ++k;
    }

    return DataFrame::create(
      _["from"] = from,
      _["to"] = to,
      _["boundary_length"] = border
    );

  } else {

    int k = 0;

    for (const auto& kv : edge_lengths) {

      uint64_t key = kv.first;

      from[k] = static_cast<int>(key >> 32);
      to[k] = static_cast<int>(key & 0xffffffff);

      ++k;
    }

    return DataFrame::create(
      _["from"] = from,
      _["to"] = to
    );
  }
}
