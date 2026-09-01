#include <Rcpp.h>
#include <nanoflann.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace nanoflann;
struct PointCloud {
  const NumericVector& x;
  const NumericVector& y;

  PointCloud(const NumericVector& x_, const NumericVector& y_)
    : x(x_), y(y_) {}

  inline size_t kdtree_get_point_count() const {
    return x.size();
  }

  inline double kdtree_get_pt(const size_t idx, int dim) const {
    return (dim == 0) ? x[idx] : y[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector cpp_idw_planar_kdtree(NumericMatrix target_xy,
                                    NumericVector obs_x,
                                    NumericVector obs_y,
                                    NumericVector obs_z,
                                    int k, double p, double rmax,
                                    bool rigorous, int cores) {

  int n_target = target_xy.nrow();
  NumericVector out(n_target);

  double rmax2 = rmax * rmax;

  // Build KD-tree
  PointCloud cloud(obs_x, obs_y);

  typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    2
  > kd_tree_t;

  kd_tree_t index(2, cloud, KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

#ifdef _OPENMP
  if (cores > 0) omp_set_num_threads(cores);
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n_target; i++) {

    double query_pt[2] = { target_xy(i, 0), target_xy(i, 1) };

    std::vector<nanoflann::ResultItem<unsigned int, double>> matches;

    nanoflann::SearchParameters params;
    params.sorted = false;
    index.radiusSearch(&query_pt[0], rmax2, matches, nanoflann::SearchParameters());
    if (k > 0 && matches.size() > (size_t)k) {
      std::nth_element(
        matches.begin(),
        matches.begin() + k,
        matches.end(),
        [](const nanoflann::ResultItem<unsigned int, double>& a,
           const nanoflann::ResultItem<unsigned int, double>& b) {
          return a.second < b.second;
        }
      );
      matches.resize(k);
    }

    double num = 0.0, den = 0.0;
    int valid_neighbors = 0;
    bool exact_match = false;

    for (size_t j = 0; j < matches.size(); j++) {
      size_t idx = matches[j].first;
      double d2 = matches[j].second;

      valid_neighbors++;

      if (d2 < 1e-12) {
        out[i] = obs_z[idx];
        exact_match = true;
        break;
      }

      double weight = 1.0 / std::pow(d2, p / 2.0);
      num += obs_z[idx] * weight;
      den += weight;
    }

    if (!exact_match) {
      if (valid_neighbors == 0 || (rigorous && valid_neighbors < k)) {
        out[i] = NA_REAL;
      } else {
        out[i] = num / den;
      }
    }
  }

  return out;
}
