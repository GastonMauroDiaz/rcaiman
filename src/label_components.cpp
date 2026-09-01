#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix label_components_cpp(LogicalMatrix bin) {

  const int nrow = bin.nrow();
  const int ncol = bin.ncol();

  IntegerMatrix labels(nrow, ncol);

  int current_label = 0;

  std::vector<int> stack;

  for (int row = 0; row < nrow; ++row) {
    for (int col = 0; col < ncol; ++col) {

      if (!bin(row, col))
        continue;

      if (labels(row, col) != 0)
        continue;

      ++current_label;

      labels(row, col) = current_label;

      stack.clear();
      stack.push_back(row * ncol + col);

      while (!stack.empty()) {

        const int idx = stack.back();
        stack.pop_back();

        const int r = idx / ncol;
        const int c = idx % ncol;

        // north
        if (r > 0 &&
            bin(r - 1, c) &&
            labels(r - 1, c) == 0) {

          labels(r - 1, c) = current_label;
          stack.push_back((r - 1) * ncol + c);
        }

        // south
        if (r < nrow - 1 &&
            bin(r + 1, c) &&
            labels(r + 1, c) == 0) {

          labels(r + 1, c) = current_label;
          stack.push_back((r + 1) * ncol + c);
        }

        // west
        if (c > 0 &&
            bin(r, c - 1) &&
            labels(r, c - 1) == 0) {

          labels(r, c - 1) = current_label;
          stack.push_back(r * ncol + (c - 1));
        }

        // east
        if (c < ncol - 1 &&
            bin(r, c + 1) &&
            labels(r, c + 1) == 0) {

          labels(r, c + 1) = current_label;
          stack.push_back(r * ncol + (c + 1));
        }
      }
    }
  }

  return labels;
}
