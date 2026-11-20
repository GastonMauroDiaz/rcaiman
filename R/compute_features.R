#' Compute object features
#'
#' @description
#' Compute object-level attributes from a labeled segmentation and reconstruct
#' raster layers from those attributes.
#'
#' @param r single-layer numeric or logical [terra::SpatRaster-class].
#'   Measurements on the grid of `seg`. Missing values are omitted from
#'   radiometric statistics; other values must be finite in included segments.
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, on the same grid as `r`. Zenith angles
#'   may span `[0, 180]`; the function does not restrict segments to the sky.
#' @param w single-layer [terra::SpatRaster-class] of area weights, typically
#'   from [area_weight_image()], on the same grid as `r`. Included weights must
#'   be finite and nonnegative, or missing. Zero and missing weights contribute
#'   no area. Relative weights suffice for means and spherical centroids.
#' @param features data frame returned by [compute_features()], or a user-created
#'   table with unique nonnegative integer `id` values and numeric or logical
#'   feature columns. Rows are matched by `id`, not by position.
#' @param field character vector. Names of columns in `features` to reconstruct
#'   as raster layers.
#'
#' @inheritParams compute_mean
#'
#' @section Functions:
#' \describe{
#'   \item{`compute_features()`}{Return a data frame with one row per segment and
#'     columns containing radiometric, geometric, and area-related attributes.}
#'   \item{`feature_image()`}{Return a numeric [terra::SpatRaster-class] whose
#'     pixel values are reconstructed from object attributes.}
#' }
#'
#' @details
#' `compute_features()` returns one row per included label, ordered by `id`.
#' Geometry describes every pixel of the segment, independently of missing
#' measurements in `r`. Radiometric means and standard deviations use `w`;
#' [compute_mean()] instead gives every nonmissing pixel equal weight.
#'
#' `feature_image()` paints each selected column over its segment. Missing
#' labels, excluded zero labels, and labels absent from `features` become `NA`.
#' Columns become layers named by `field`; unused table rows are ignored.
#'
#' Features returned by `compute_features()`:
#'
#' `id` is the original segment label.
#'
#' \describe{
#'   \item{`pixel_count`}{Number of pixels belonging to the segment.}
#'
#'   \item{`valid_count`}{Number of non-missing pixels.}
#'
#'   \item{`perimeter`}{
#'       Segment perimeter estimated from contour pixels (4 connectivity) and
#'       their associated area weights.
#'       Each contour pixel contributes a length proportional to the diameter
#'       of a circle with equivalent area:
#'       \deqn{
#'         d_i = 2 \sqrt{\frac{w_i}{\pi}}
#'       }
#'       where \eqn{w_i} is the area weight associated with pixel \eqn{i}. This
#'       quantity is a shape descriptor. It should not be interpreted as an
#'       exact geodesic perimeter on the hemisphere.
#'      }
#'
#'   \item{`area`}{
#'      Segment area computed as:
#'         \deqn{
#'           A = \sum_i w_i
#'         }
#'      where \eqn{w_i} is the area weight associated with pixel \eqn{i}. Equal
#'      weighted areas therefore represent approximate equal portions of
#'      angular space regardless of image position. Its units follow `w`:
#'      relative weights give relative area, not steradians.
#'     }
#'
#'   \item{`shape_complexity`}{Dimensionless shape descriptor computed as
#'     `perimeter^2 / area`. Larger values indicate more irregular or
#'     elongated shapes.}
#'
#'   \item{`effective_pixel_area`}{Effective object size computed as
#'     `pixel_count / shape_complexity`. Large values indicate compact regions
#'     represented by many pixels.}
#'
#'   \item{`effective_sample_size`}{
#'     Effective sample size induced by area weighting:
#'     \deqn{
#'       n_{\mathrm{eff}}
#'       =
#'       \frac{
#'         (\sum_i w_i)^2
#'       }{
#'         \sum_i w_i^2
#'       }
#'     }
#'     Used internally in the correction factor for `sd`. With equal positive
#'     weights, it equals `valid_count`; it decreases as the weighted mean
#'     becomes dominated by fewer pixels. Weighting accounts for the area
#'     represented by each pixel, but a larger weight does not supply additional
#'     measurements. A low value therefore helps identify means that depend
#'     heavily on a few pixel values, even when the segment contains many pixels.
#'    }
#'
#'   \item{`mean`}{
#'     Area-weighted mean pixel value:
#'     \deqn{
#'       \bar{x}_w =
#'       \frac{\sum_i w_i x_i}
#'            {\sum_i w_i}
#'     }
#'   }
#'
#'   \item{`sd`}{
#'     Area-weighted sample standard deviation:
#'     \deqn{
#'       s_w =
#'       \sqrt{
#'         \frac{
#'           \sum_i w_i
#'           (x_i-\bar{x}_w)^2
#'         }{
#'           \sum_i w_i
#'         }
#'         \cdot
#'         \frac{
#'           n_{\mathrm{eff}}
#'         }{
#'           n_{\mathrm{eff}} - 1
#'         }
#'       }
#'     }
#'     Here \eqn{n_{\mathrm{eff}}} is `effective_sample_size`. The correction
#'     factor is applied automatically; it requires no separate user input.
#'     This estimator reduces to the ordinary sample standard deviation when
#'     all pixels have equal weight.
#'   }
#'
#'   \item{`min`}{Minimum pixel value.}
#'
#'   \item{`max`}{Maximum pixel value.}
#'
#'   \item{`row`, `col`}{Coordinates in raster
#'     space of the pixel with the maximum value.}
#'
#'   \item{`zenith`, `azimuth`}{External direction of the maximum in radians,
#'     as in [add_directions()]. Ties select the first pixel in raster cell order.}
#'
#'   \item{`spherical_centroid_zenith`, `spherical_centroid_azimuth`}{
#'     Direction of the area-weighted centroid on the unit sphere, in radians.
#'     The three-dimensional unit vectors are averaged using `w`, and the mean
#'     vector is normalized to unit length. Area weights account for differences
#'     in the solid angle represented by pixels under the image projection;
#'     expressing directions as unit vectors alone does not correct this effect.
#'     Relative area weights suffice: multiplying all weights by the same positive
#'     constant leaves the direction unchanged. Missing values in `r` do not
#'     affect this geometric descriptor. Zero or missing weights do not contribute;
#'     a missing direction with positive weight makes the centroid undefined.
#'     Both angles are `NA` when total weight is zero or the mean vector has
#'     numerically zero length.}
#'
#'   \item{`planar_centroid_row`, `planar_centroid_col`}{
#'     Geometric centroid of the segment in the image plane. The centroid is
#'     obtained by simple averaging.
#'   }
#' }
#' @return See *Functions*.
#'
#' @seealso [compute_mean()]
#'
#' @name compute_features
#' @rdname compute_features
#' @aliases feature_image
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' w <- area_weight_image(g, geometry)
#' seg <- sky_grid_segmentation(d, 15)
#' # Illustrative measurements; supply an image band in a real analysis.
#' r <- cos(d$zenith * pi / 180)
#' features <- compute_features(r, seg, d, w)
#' head(features)
#' plot(feature_image(seg, features, c("mean", "area")))
compute_features <- function(r, seg, d, w, ignore_label_0 = TRUE) {
  .assert_single_layer(r)
  .assert_same_geom(r, seg)
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  .assert_same_geom(r, z)
  .assert_single_layer(w)
  .assert_same_geom(r, w)
  labels <- .feature_labels(seg, ignore_label_0)
  included <- !is.na(labels$index)
  values <- terra::values(r, mat = FALSE)
  weights <- terra::values(w, mat = FALSE)
  .check_vector(as.numeric(values[included]), "numeric", allow_na = TRUE, name = "r")
  .check_vector(weights[included], "numeric", allow_na = TRUE,
                sign = "nonnegative", name = "w")
  f <- compute_features_cpp(
    values = values,
    labels = labels$index,
    zenith = terra::values(z, mat = FALSE) * pi / 180,
    azimuth = terra::values(a, mat = FALSE) * pi / 180,
    weights = weights,
    nrow = nrow(r), ncol = ncol(r), ignore_label_0 = FALSE
  )
  f$id <- labels$ids[f$id]
  f$shape_complexity <- ifelse(f$area > 0, f$perimeter^2 / f$area, NA_real_)
  f$effective_pixel_area <- f$pixel_count / f$shape_complexity
  f
}

#' @rdname compute_features
#' @export
feature_image <- function(seg, features, field, ignore_label_0 = TRUE) {
  labels <- .feature_labels(seg, ignore_label_0)
  .check_vector(field, "character")
  if (!length(field) || anyDuplicated(field)) {
    stop("`field` must contain distinct column names.", call. = FALSE)
  }
  .check_columns(features, required_columns = c("id", field))
  .check_vector(features$id, "integerish", sign = "nonnegative", name = "features$id")
  if (anyDuplicated(features$id)) stop("`features$id` must be unique.", call. = FALSE)
  if (!all(vapply(features[field], function(x) is.numeric(x) || is.logical(x), logical(1)))) {
    stop("Feature columns must be numeric or logical.", call. = FALSE)
  }
  rows <- match(labels$ids[labels$index], features$id)
  values <- as.matrix(features[, field, drop = FALSE])[rows, , drop = FALSE]
  out <- terra::rast(seg, nlyrs = length(field))
  terra::values(out) <- values
  names(out) <- field
  out
}

# Compact labels for native routines while preserving the original IDs in R.
.feature_labels <- function(seg, ignore_label_0) {
  .assert_single_layer(seg)
  .check_vector(ignore_label_0, "logical", 1)
  values <- .segment_values(seg)
  ids <- sort(unique(values[!is.na(values) & (!ignore_label_0 | values != 0)]))
  list(ids = ids, index = match(values, ids))
}
