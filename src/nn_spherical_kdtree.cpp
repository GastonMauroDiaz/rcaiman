// file: nn_spherical_kdtree.cpp

#include <Rcpp.h>
#include <nanoflann.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace nanoflann;


// ============================================================================
// Point cloud
// ============================================================================

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


// ============================================================================
// Main
// ============================================================================

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List nn_spherical_cpp(NumericVector sp_x,
                      NumericVector sp_y,
                      NumericVector sp_z,
                      NumericVector qp_x,
                      NumericVector qp_y,
                      NumericVector qp_z,
                      int k = 1,
                      Nullable<double> search_radius = R_NilValue,
                      int cores = 1) {

  int n_sp = sp_x.size();
  int n_qp = qp_x.size();

  if (sp_y.size() != n_sp || sp_z.size() != n_sp) {
    stop("Sampling point vectors must have same length.");
  }

  if (qp_y.size() != n_qp || qp_z.size() != n_qp) {
    stop("Query point vectors must have same length.");
  }

  if (k < 1) {
    stop("`k` must be positive.");
  }

  if (k > n_sp) {
    k = n_sp;
  }


  // ==========================================================================
  // Build KD-tree
  // ==========================================================================

  PointCloud3D cloud(sp_x, sp_y, sp_z);

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


  // ==========================================================================
  // Radius conversion
  // ==========================================================================

  bool use_radius = search_radius.isNotNull();

  double radius_chord2 = NA_REAL;

  if (use_radius) {

    double theta = as<double>(search_radius);

    // chord^2 = 2 - 2*cos(theta)
    radius_chord2 = 2.0 - 2.0 * std::cos(theta);
  }


  // ==========================================================================
  // Output
  // ==========================================================================

  List output(n_qp);


#ifdef _OPENMP
  if (cores > 0) omp_set_num_threads(cores);
#endif


  // ==========================================================================
  // Query loop
  // ==========================================================================

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n_qp; i++) {

    double query_pt[3] = {
      qp_x[i],
          qp_y[i],
              qp_z[i]
    };


    // ========================================================================
    // KNN search
    // ========================================================================

    std::vector<size_t> ret_index(k);
    std::vector<double> out_dist_sqr(k);

    KNNResultSet<double> resultSet(k);

    resultSet.init(
      ret_index.data(),
      out_dist_sqr.data()
    );

    index.findNeighbors(
      resultSet,
      &query_pt[0],
               nanoflann::SearchParameters()
    );


    // ========================================================================
    // Convert distances
    // ========================================================================

    std::vector<int> out_index;
    std::vector<double> out_distance;

    out_index.reserve(k);
    out_distance.reserve(k);

    for (int j = 0; j < k; j++) {

      double d2 = out_dist_sqr[j];

      if (use_radius && d2 > radius_chord2) {
        continue;
      }

      size_t idx = ret_index[j];

      double dot =
        sp_x[idx] * query_pt[0] +
        sp_y[idx] * query_pt[1] +
        sp_z[idx] * query_pt[2];

      if (dot > 1.0) {
        dot = 1.0;
      } else if (dot < -1.0) {
        dot = -1.0;
      }

      double angular_distance = std::acos(dot);

      out_index.push_back(idx + 1);
      out_distance.push_back(angular_distance);
    }

#ifdef _OPENMP
#pragma omp critical
#endif
{
  output[i] = List::create(
    Named("index") = out_index,
    Named("distance") = out_distance
  );
}
  }

  return output;
}
