#include <Rcpp.h>
#include <unordered_map>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix feature_image_cpp(
    IntegerVector labels,
    IntegerVector ids,
    NumericMatrix feature,
    bool ignore_label_0 = true
) {

  const R_xlen_t n = labels.size();

  if (ids.size() != feature.nrow()) {
    stop("`ids` length must equal the number of rows in `feature`.");
  }

  const int n_features = feature.ncol();

  // label -> row index in feature matrix
  std::unordered_map<int, int> lookup;
  lookup.reserve(ids.size());

  for (R_xlen_t i = 0; i < ids.size(); ++i) {

    int id = ids[i];

    if (IntegerVector::is_na(id))
      continue;

    lookup[id] = static_cast<int>(i);
  }

  NumericMatrix out(n, n_features);

  std::fill(
    out.begin(),
    out.end(),
    NA_REAL
  );

  for (R_xlen_t i = 0; i < n; ++i) {

    int label = labels[i];

    if (IntegerVector::is_na(label))
      continue;

    if (ignore_label_0 && label == 0)
      continue;

    auto it = lookup.find(label);

    if (it == lookup.end())
      continue;

    int row_feature = it->second;

    for (int j = 0; j < n_features; ++j) {
      out(i, j) = feature(row_feature, j);
    }
  }

  return out;
}
