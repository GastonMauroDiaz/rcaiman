#include <Rcpp.h>
#include <unordered_map>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List compute_mean_cpp(NumericVector values,
                            IntegerVector labels,
                            bool ignore_label_0 = true) {

  std::unordered_map<int, double> sums;
  std::unordered_map<int, int> counts;

  const R_xlen_t n = values.size();

  for (R_xlen_t i = 0; i < n; ++i) {

    if (NumericVector::is_na(values[i]) ||
        IntegerVector::is_na(labels[i])) {
      continue;
    }

    int label = static_cast<int>(labels[i]);

    if (ignore_label_0 && label == 0) {
      continue;
    }

    sums[label] += values[i];
    counts[label] += 1;
  }

  const std::size_t n_labels = sums.size();

  NumericVector means(n_labels);
  IntegerVector ids(n_labels);

  std::size_t k = 0;

  for (const auto& pair : sums) {
    int label = pair.first;

    means[k] = pair.second / counts[label];
    ids[k] = label;

    ++k;
  }

  return DataFrame::create(
    _["id"] = ids,
    _["mean"] = means
  );
}
