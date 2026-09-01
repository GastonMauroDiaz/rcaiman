// file: src/thr_nbhu.cpp

#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List thr_nbhu_cpp(
    NumericMatrix img,
    int min_thresh = 0,
    int max_thresh = 254
) {

  const int nrows = img.nrow();
  const int ncols = img.ncol();

  const int n_thresholds =
    max_thresh - min_thresh + 1;

  NumericVector scores(n_thresholds);
  IntegerVector ntv(n_thresholds);
  IntegerVector thresholds(n_thresholds);

  int idx = 0;

  for (int t = min_thresh; t <= max_thresh; ++t) {

    double sum_diff = 0.0;
    int n_transitions = 0;

    for (int r = 0; r < (nrows - 1); ++r) {

      for (int c = 0; c < (ncols - 1); ++c) {

        const int p11 =
          static_cast<int>(img(r, c));

        const int p12 =
          static_cast<int>(img(r, c + 1));

        const int p21 =
          static_cast<int>(img(r + 1, c));

        const int p22 =
          static_cast<int>(img(r + 1, c + 1));

        const bool b11 = (p11 > t);
        const bool b12 = (p12 > t);
        const bool b21 = (p21 > t);
        const bool b22 = (p22 > t);

        // 11 vs 12
        if (b11 != b12) {
          sum_diff += std::abs(p11 - p12);
          ++n_transitions;
        }

        // 11 vs 21
        if (b11 != b21) {
          sum_diff += std::abs(p11 - p21);
          ++n_transitions;
        }

        // 11 vs 22
        if (b11 != b22) {
          sum_diff += std::abs(p11 - p22);
          ++n_transitions;
        }

        // 12 vs 21
        if (b12 != b21) {
          sum_diff += std::abs(p12 - p21);
          ++n_transitions;
        }
      }
    }

    double score = NA_REAL;

    if (n_transitions > 0) {

      score =
        sum_diff /
          static_cast<double>(n_transitions);
    }

    thresholds[idx] = t;
    scores[idx] = score;
    ntv[idx] = n_transitions;

    ++idx;
  }

  return List::create(
    _["thr"] = thresholds,
    _["score"] = scores,
    _["n"] = ntv
  );
}
