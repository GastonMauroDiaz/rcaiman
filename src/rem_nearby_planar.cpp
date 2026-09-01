// rem_nearby_planar.cpp

#include <Rcpp.h>
#include <nanoflann.hpp>

using namespace Rcpp;
using namespace nanoflann;

struct PointCloud2D {
  const NumericVector& x;
  const NumericVector& y;

  PointCloud2D(const NumericVector& x_,
               const NumericVector& y_)
    : x(x_), y(y_) {}

  inline size_t kdtree_get_point_count() const {
    return x.size();
  }

  inline double kdtree_get_pt(const size_t idx, int dim) const {
    return (dim == 0) ? x[idx] : y[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const {
    return false;
  }
};

// [[Rcpp::export]]
LogicalVector cpp_rem_nearby_planar(
    NumericVector x,
    NumericVector y,
    double min_dist
) {

  size_t n = x.size();

  LogicalVector keep(n, false);
  std::vector<bool> discarded(n, false);

  PointCloud2D cloud(x, y);

  typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud2D>,
    PointCloud2D,
    2
  > kd_tree_t;

  kd_tree_t index(
      2,
      cloud,
      KDTreeSingleIndexAdaptorParams(10)
  );

  index.buildIndex();

  double radius2 = min_dist * min_dist;

  for (size_t i = 0; i < n; ++i) {

    if (discarded[i]) continue;

    keep[i] = true;

    double query_pt[2] = { x[i], y[i] };

    std::vector<nanoflann::ResultItem<unsigned int, double>> matches;

    nanoflann::SearchParameters params;
    params.sorted = false;

    index.radiusSearch(
      &query_pt[0],
               radius2,
               matches,
               params
    );

    for (size_t j = 0; j < matches.size(); ++j) {
      discarded[matches[j].first] = true;
    }

    discarded[i] = false;
  }

  return keep;
}
