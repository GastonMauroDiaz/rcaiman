#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

struct FeatureAccumulator {

  uint64_t pixel_count = 0;
  uint64_t valid_count = 0;

  double area = 0.0;

  // double sum_value = 0.0;
  // double sum_value_sq = 0.0;

  double sum_weighted_value = 0.0;
  double sum_weighted_value_sq = 0.0;

  double valid_weight = 0.0;
  double valid_weight_sq = 0.0;

  double min_value = R_PosInf;
  double max_value = R_NegInf;

  double max_row = NA_REAL;
  double max_col = NA_REAL;

  double max_zenith = NA_REAL;
  double max_azimuth = NA_REAL;

  double sum_row = 0.0;
  double sum_col = 0.0;

  double sum_weighted_row = 0.0;
  double sum_weighted_col = 0.0;

  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_z = 0.0;

  double sum_xx = 0.0;
  double sum_yy = 0.0;
  double sum_zz = 0.0;

  double sum_xy = 0.0;
  double sum_xz = 0.0;
  double sum_yz = 0.0;

  double perimeter = 0.0;
};

// [[Rcpp::export]]
DataFrame compute_features_cpp(
    NumericVector values,
    IntegerVector labels,
    NumericVector zenith,
    NumericVector azimuth,
    Nullable<NumericVector> weights = R_NilValue,
    int nrow = 0,
    int ncol = 0,
    bool ignore_label_0 = true
) {

  const R_xlen_t n = labels.size();

  if (values.size() != n) {
    stop("`values` and `labels` must have the same length.");
  }

  if (zenith.size() != n) {
    stop("`zenith` must have the same length as `labels`.");
  }

  if (azimuth.size() != n) {
    stop("`azimuth` must have the same length as `labels`.");
  }

  bool use_weights = weights.isNotNull();

  NumericVector w;

  if (use_weights) {
    w = NumericVector(weights);
    if (w.size() != n) {
      stop("`weights` must have the same length as `labels`.");
    }
  }

  if ((nrow * ncol) != n) {
    stop("`nrow * ncol` must equal the number of pixels.");
  }

  // -------------------------------------------------------------------------
  // Extract unique labels
  // -------------------------------------------------------------------------

  std::vector<int> ids;
  ids.reserve(n / 10);

  for (R_xlen_t i = 0; i < n; ++i) {

    int label = labels[i];

    if (IntegerVector::is_na(label))
      continue;

    if (ignore_label_0 && label == 0)
      continue;

    ids.push_back(label);
  }

  std::sort(ids.begin(), ids.end());

  ids.erase(
    std::unique(ids.begin(), ids.end()),
    ids.end()
  );

  const size_t n_segments = ids.size();

  std::unordered_map<int, size_t> id_to_index;
  id_to_index.reserve(n_segments);

  for (size_t i = 0; i < n_segments; ++i) {
    id_to_index[ids[i]] = i;
  }

  std::vector<FeatureAccumulator> acc(n_segments);

  // -------------------------------------------------------------------------
  // Main loop
  // -------------------------------------------------------------------------

  for (R_xlen_t i = 0; i < n; ++i) {

    int label = labels[i];

    if (IntegerVector::is_na(label))
      continue;

    if (ignore_label_0 && label == 0)
      continue;

    size_t idx = id_to_index[label];

    FeatureAccumulator& a = acc[idx];

    double weight = 1.0;

    if (use_weights) {

      weight = w[i];

      if (NumericVector::is_na(weight))
        weight = 0.0;
    }

    // row/col coordinates (1-based)

    double row = static_cast<double>(i / ncol) + 1.0;
    double col = static_cast<double>(i % ncol) + 1.0;

    double theta = zenith[i];
    double phi = azimuth[i];

    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);

    a.pixel_count++;

    // START perimeter
    bool is_boundary = false;

    int r = i / ncol;
    int c = i % ncol;

    // north
    if (r == 0) {

      is_boundary = true;

    } else {

      int nb = labels[i - ncol];

      if (IntegerVector::is_na(nb) || nb != label)
        is_boundary = true;
    }

    // south
    if (!is_boundary) {

      if (r == (nrow - 1)) {

        is_boundary = true;

      } else {

        int nb = labels[i + ncol];

        if (IntegerVector::is_na(nb) || nb != label)
          is_boundary = true;
      }
    }

    // west
    if (!is_boundary) {

      if (c == 0) {

        is_boundary = true;

      } else {

        int nb = labels[i - 1];

        if (IntegerVector::is_na(nb) || nb != label)
          is_boundary = true;
      }
    }

    // east
    if (!is_boundary) {

      if (c == (ncol - 1)) {

        is_boundary = true;

      } else {

        int nb = labels[i + 1];

        if (IntegerVector::is_na(nb) || nb != label)
          is_boundary = true;
      }
    }

    if (is_boundary) {

      a.perimeter +=
        2.0 * std::sqrt(weight / M_PI);
    }
    // END perimeter

    a.area += weight;

    a.sum_row += row;
    a.sum_col += col;

    a.sum_weighted_row += row * weight;
    a.sum_weighted_col += col * weight;

    // Unit vectors describe directions; weights account for pixel solid angle.
    // A zero-area pixel contributes nothing, even if its direction is missing.
    if (weight > 0.0) {
      a.sum_x += weight * x;
      a.sum_y += weight * y;
      a.sum_z += weight * z;
    }

    a.sum_xx += x * x;
    a.sum_yy += y * y;
    a.sum_zz += z * z;

    a.sum_xy += x * y;
    a.sum_xz += x * z;
    a.sum_yz += y * z;

    double value = values[i];

    if (!NumericVector::is_na(value)) {

      a.valid_count++;

      // a.sum_value += value;
      // a.sum_value_sq += value * value;

      a.sum_weighted_value += weight * value;
      a.sum_weighted_value_sq += weight * value * value;

      a.valid_weight += weight;
      a.valid_weight_sq += weight * weight;

      if (value < a.min_value)
        a.min_value = value;

      if (value > a.max_value) {

        a.max_value = value;

        a.max_row = row;
        a.max_col = col;

        a.max_zenith = theta;
        a.max_azimuth = phi;
      }
    }
  }

  // -------------------------------------------------------------------------
  // Build output
  // -------------------------------------------------------------------------

  IntegerVector out_id(n_segments);

  NumericVector pixel_count(n_segments);
  NumericVector valid_count(n_segments);

  NumericVector area(n_segments);

  NumericVector mean(n_segments);
  NumericVector sd(n_segments);

  NumericVector effective_sample_size(n_segments);

  NumericVector min_value(n_segments);
  NumericVector max_value(n_segments);

  NumericVector max_row(n_segments);
  NumericVector max_col(n_segments);

  NumericVector max_zenith(n_segments);
  NumericVector max_azimuth(n_segments);

  NumericVector centroid_row(n_segments);
  NumericVector centroid_col(n_segments);

  NumericVector weighted_centroid_row(n_segments);
  NumericVector weighted_centroid_col(n_segments);

  NumericVector centroid_zenith(n_segments);
  NumericVector centroid_azimuth(n_segments);

  NumericVector perimeter(n_segments);

  for (size_t i = 0; i < n_segments; ++i) {

    const FeatureAccumulator& a = acc[i];

    out_id[i] = ids[i];

    pixel_count[i] = static_cast<double>(a.pixel_count);
    valid_count[i] = static_cast<double>(a.valid_count);

    area[i] = a.area;

    max_row[i] = a.max_row;
    max_col[i] = a.max_col;

    max_zenith[i] = a.max_zenith;
    max_azimuth[i] = a.max_azimuth;

    centroid_row[i] =
      a.pixel_count > 0 ?
      a.sum_row / a.pixel_count :
      NA_REAL;

    centroid_col[i] =
      a.pixel_count > 0 ?
      a.sum_col / a.pixel_count :
      NA_REAL;

    weighted_centroid_row[i] =
      a.area > 0.0 ?
      a.sum_weighted_row / a.area :
      NA_REAL;

    weighted_centroid_col[i] =
      a.area > 0.0 ?
      a.sum_weighted_col / a.area :
      NA_REAL;

    perimeter[i] = a.perimeter;

    if (a.area > 0.0) {

      double mx = a.sum_x / a.area;
      double my = a.sum_y / a.area;
      double mz = a.sum_z / a.area;

      double norm =
        std::sqrt(
          mx * mx +
            my * my +
            mz * mz
        );

      if (norm > 64.0 * std::numeric_limits<double>::epsilon()) {

        mx /= norm;
        my /= norm;
        mz /= norm;

        centroid_zenith[i] =
          std::acos(std::max(-1.0, std::min(1.0, mz)));

        double az =
          std::atan2(my, mx);

        if (az < 0.0)
          az += 2.0 * M_PI;

        centroid_azimuth[i] = az;

      } else {

        centroid_zenith[i] = NA_REAL;
        centroid_azimuth[i] = NA_REAL;
      }

    } else {

      centroid_zenith[i] = NA_REAL;
      centroid_azimuth[i] = NA_REAL;
    }

    min_value[i] = a.valid_count > 0 ? a.min_value : NA_REAL;
    max_value[i] = a.valid_count > 0 ? a.max_value : NA_REAL;

    if (a.valid_count > 0 &&
        a.valid_weight > 0.0) {

      double m =
        a.sum_weighted_value /
          a.valid_weight;

      mean[i] = m;

      min_value[i] = a.min_value;
      max_value[i] = a.max_value;

      double n_eff =
        (a.valid_weight * a.valid_weight) /
          a.valid_weight_sq;

      effective_sample_size[i] = n_eff;

      if (n_eff > 1.0) {

        double variance_pop =
          a.sum_weighted_value_sq /
            a.valid_weight -
              m * m;

        if (variance_pop < 0.0)
          variance_pop = 0.0;

        double variance =
          variance_pop *
          n_eff /
            (n_eff - 1.0);

        sd[i] = std::sqrt(variance);

      } else {

        sd[i] = NA_REAL;
      }

    } else {

      mean[i] = NA_REAL;
      sd[i] = NA_REAL;

      effective_sample_size[i] = NA_REAL;


    }
  }

  return DataFrame::create(
    _["id"] = out_id,
    _["pixel_count"] = pixel_count,
    _["valid_count"] = valid_count,
    _["effective_sample_size"] = effective_sample_size,
    _["perimeter"] = perimeter,
    _["area"] = area,
    _["mean"] = mean,
    _["sd"] = sd,
    _["min"] = min_value,
    _["max"] = max_value,
    _["row"] = max_row,
    _["col"] = max_col,
    _["zenith"] = max_zenith,
    _["azimuth"] = max_azimuth,
    _["planar_centroid_row"] = centroid_row,
    _["planar_centroid_col"] = centroid_col,
    // _["weighted_centroid_row"] = weighted_centroid_row,
    // _["weighted_centroid_col"] = weighted_centroid_col,
    _["spherical_centroid_zenith"] = centroid_zenith,
    _["spherical_centroid_azimuth"] = centroid_azimuth,
    _["stringsAsFactors"] = false
  );
}
