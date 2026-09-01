// rem_nearby_spherical.cpp

#include <Rcpp.h>
#include <nanoflann.hpp>

using namespace Rcpp;
using namespace nanoflann;

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

// [[Rcpp::export]]
LogicalVector cpp_rem_nearby_spherical(
    NumericVector x,
    NumericVector y,
    NumericVector z,
    double min_angle_rad
) {

  size_t n = x.size();

  LogicalVector keep(n, false);
  std::vector<bool> discarded(n, false);

  PointCloud3D cloud(x, y, z);

  typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud3D>,
    PointCloud3D,
    3
  > kd_tree_t;

  kd_tree_t index(
      3,
      cloud,
      KDTreeSingleIndexAdaptorParams(10)
  );

  index.buildIndex();

  double chord = 2.0 * std::sin(min_angle_rad / 2.0);
  double radius2 = chord * chord;

  for (size_t i = 0; i < n; ++i) {

    if (discarded[i]) continue;

    keep[i] = true;

    double query_pt[3] = {
      x[i],
       y[i],
        z[i]
    };

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
