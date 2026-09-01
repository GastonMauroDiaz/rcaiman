#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector robust_cv(
    NumericVector x,
    size_t ni,
    size_t nw
) {

  NumericVector out(ni);

  size_t start = 0;

  for (size_t i=0; i<ni; i++) {

    size_t end = start + nw;

    double vals[9];
    int n = 0;

    for (size_t j=start; j<end; j++) {
      if (!NumericVector::is_na(x[j])) {
        vals[n++] = x[j];
      }
    }

    if (n == 0) {
      out[i] = NA_REAL;
      start = end;
      continue;
    }

    std::nth_element(
      vals,
      vals + n/2,
      vals + n
    );

    double med = vals[n/2];

    if (n % 2 == 0) {
      double upper = med;

      std::nth_element(
        vals,
        vals + n/2 - 1,
        vals + n
      );

      med = (upper + vals[n/2 -1]) / 2.0;
    }

    if (med == 0) {
      out[i] = NA_REAL;
      start = end;
      continue;
    }

    for (int k=0; k<n; k++) {
      vals[k] = std::abs(vals[k] - med);
    }

    std::nth_element(
      vals,
      vals + n/2,
      vals + n
    );

    double mad = vals[n/2];

    if (n % 2 == 0) {
      double upper = mad;

      std::nth_element(
        vals,
        vals + n/2 -1,
        vals + n
      );

      mad = (upper + vals[n/2 -1]) / 2.0;
    }

    out[i] = (1.4826 * mad) / med;

    start = end;
  }

  return out;
}
