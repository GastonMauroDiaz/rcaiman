#include <Rcpp.h>
#include <vector>
#include <queue>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix distance_to_black_cpp(LogicalMatrix bin) {

  const int nrow = bin.nrow();
  const int ncol = bin.ncol();

  IntegerMatrix dist(nrow, ncol);

  std::queue<int> q;

  bool has_black = false;

  // Initialize
  for (int row = 0; row < nrow; ++row) {
    for (int col = 0; col < ncol; ++col) {

      if (bin(row, col)) {

        dist(row, col) = -1;

      } else {

        dist(row, col) = 0;

        q.push(row * ncol + col);

        has_black = true;
      }
    }
  }

  if (!has_black) {

    warning(
      "`bin` does not contain any FALSE pixels. "
      "Distance to black is undefined. Returning NA."
    );

    IntegerMatrix out(nrow, ncol);

    std::fill(out.begin(), out.end(), NA_INTEGER);

    return out;
  }

  while (!q.empty()) {

    const int idx = q.front();
    q.pop();

    const int r = idx / ncol;
    const int c = idx % ncol;

    const int next_dist = dist(r, c) + 1;

    // NW
    if (r > 0 && c > 0 &&
        dist(r - 1, c - 1) == -1) {

      dist(r - 1, c - 1) = next_dist;
      q.push((r - 1) * ncol + (c - 1));
    }

    // N
    if (r > 0 &&
        dist(r - 1, c) == -1) {

      dist(r - 1, c) = next_dist;
      q.push((r - 1) * ncol + c);
    }

    // NE
    if (r > 0 && c < ncol - 1 &&
        dist(r - 1, c + 1) == -1) {

      dist(r - 1, c + 1) = next_dist;
      q.push((r - 1) * ncol + (c + 1));
    }

    // W
    if (c > 0 &&
        dist(r, c - 1) == -1) {

      dist(r, c - 1) = next_dist;
      q.push(r * ncol + (c - 1));
    }

    // E
    if (c < ncol - 1 &&
        dist(r, c + 1) == -1) {

      dist(r, c + 1) = next_dist;
      q.push(r * ncol + (c + 1));
    }

    // SW
    if (r < nrow - 1 && c > 0 &&
        dist(r + 1, c - 1) == -1) {

      dist(r + 1, c - 1) = next_dist;
      q.push((r + 1) * ncol + (c - 1));
    }

    // S
    if (r < nrow - 1 &&
        dist(r + 1, c) == -1) {

      dist(r + 1, c) = next_dist;
      q.push((r + 1) * ncol + c);
    }

    // SE
    if (r < nrow - 1 && c < ncol - 1 &&
        dist(r + 1, c + 1) == -1) {

      dist(r + 1, c + 1) = next_dist;
      q.push((r + 1) * ncol + (c + 1));
    }
  }

  return dist;
}
