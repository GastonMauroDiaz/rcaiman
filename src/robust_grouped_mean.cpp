// file: robust_grouped_mean.cpp

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace Rcpp;

struct Record {
  int cell;
  double value;
};

double compute_median(std::vector<double>& x) {

  size_t n = x.size();

  if (n == 0) {
    return NA_REAL;
  }

  size_t mid = n / 2;

  std::nth_element(x.begin(), x.begin() + mid, x.end());

  double med = x[mid];

  if (n % 2 == 0) {
    std::nth_element(x.begin(), x.begin() + mid - 1, x.end());
    med = (med + x[mid - 1]) / 2.0;
  }

  return med;
}

double compute_mad(const std::vector<double>& x, double center) {

  size_t n = x.size();

  if (n == 0 || NumericVector::is_na(center)) {
    return NA_REAL;
  }

  std::vector<double> dev(n);

  for (size_t i = 0; i < n; ++i) {
    dev[i] = std::abs(x[i] - center);
  }

  double mad = compute_median(dev);

  if (NumericVector::is_na(mad)) {
    return NA_REAL;
  }

  return mad * 1.4826;
}

// [[Rcpp::export]]
DataFrame robust_grouped_mean_cpp(IntegerVector cell,
                                  NumericVector value,
                                  double laxity = 2.5) {

  int n = cell.size();

  if (value.size() != n) {
    stop("`cell` and `value` must have same length.");
  }

  std::vector<Record> records(n);

  for (int i = 0; i < n; ++i) {
    records[i] = {cell[i], value[i]};
  }

  std::sort(records.begin(), records.end(),
            [](const Record& a, const Record& b) {
              return a.cell < b.cell;
            });

  std::vector<int> out_cell;
  std::vector<double> out_value;
  std::vector<int> out_n;

  int start = 0;

  while (start < n) {

    int current_cell = records[start].cell;

    int end = start;

    while (end < n && records[end].cell == current_cell) {
      ++end;
    }

    std::vector<double> vals;
    vals.reserve(end - start);

    for (int i = start; i < end; ++i) {

      double v = records[i].value;

      if (!NumericVector::is_na(v)) {
        vals.push_back(v);
      }
    }

    int count = vals.size();

    if (count > 0) {

      std::vector<double> vals_copy = vals;

      double med = compute_median(vals_copy);

      double mad = compute_mad(vals, med);

      double sum = 0.0;
      int kept = 0;

      if (mad == 0 || NumericVector::is_na(mad)) {

        for (double v : vals) {
          sum += v;
          ++kept;
        }

      } else {

        for (double v : vals) {

          double dev = std::abs(v - med) / mad;

          if (dev <= laxity || NumericVector::is_na(dev)) {
            sum += v;
            ++kept;
          }
        }
      }

      if (kept > 0) {

        out_cell.push_back(current_cell);
        out_value.push_back(sum / kept);
        out_n.push_back(count);
      }
    }

    start = end;
  }

  return DataFrame::create(
    Named("cell") = out_cell,
    Named("value") = out_value,
    Named("n") = out_n
  );
}
