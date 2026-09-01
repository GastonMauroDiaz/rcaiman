#' Compute object features
#'
#' @description
#' Compute object-level attributes from a labeled segmentation and reconstruct
#' raster layers from those attributes.
#'
#' @param features data frame returned by [compute_features()].
#' @param field character vector. Names of columns in `features` to reconstruct
#'   as raster layers.
#'
#' @inheritParams add_directions
#' @inheritParams merge_segments
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
#' Features returned by `compute_features()`:
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
#'      hemispherical space regardless of image position.
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
#'     This quantity represents the number of equally weighted observations
#'     providing the same statistical information as the weighted sample.
#'     It equals `valid_count` when all pixels have identical weights and
#'     decreases as weight inequality increases.
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
#'   \item{`zen`, `azi`}{Direction of the maximum value as in [add_directions()].}
#'
#'   \item{`spherical_centroid_zenith`, `spherical_centroid_azimuth`}{
#'     Geometric centroid of the segment on the unit hemisphere. The centroid is
#'     obtained by averaging the three-dimensional unit vectors associated with all
#'     segment pixels and converting the resulting mean direction back to zenith and
#'     azimuth angles (radian).}
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
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' print(compute_mean(r$Cyan, g, return = "vector"))
#' plot(compute_mean(r$Cyan, g, return = "raster"))
#' }
compute_features <- function(r, seg, z, a, w, ignore_label_0 = TRUE) {

  .assert_single_layer(r, allow_na = TRUE)
  .assert_single_layer(seg)
  .assert_same_geom(r, seg)
  # .assert_zenith_image(z)
  .assert_single_layer(z, squared = FALSE)
  .assert_raster_range(z, 0, 91, open_upper = FALSE)
  .assert_same_geom(z, r)
  # .assert_azimuth_image(a)
  .assert_single_layer(a, squared = FALSE)
  .assert_raster_range(a, 0, 360, open_upper = TRUE)
  .assert_same_geom(z, a)
  .assert_single_layer(w, allow_na = TRUE)
  .assert_same_geom(z, w)
  .check_vector(ignore_label_0, "logical", 1)

  f <- compute_features_cpp(
    values = terra::values(r, mat = FALSE),
    labels = terra::as.matrix(seg, mat = FALSE),
    zenith = terra::values(z, mat = FALSE) * pi / 180,
    azimuth = terra::values(a, mat = FALSE) * pi / 180,
    weights = terra::as.matrix(w, mat = FALSE),
    nrow = nrow(r),
    ncol = ncol(r),
    ignore_label_0 = ignore_label_0
  )

  f$shape_complexity <- f$perimeter^2 / f$area
  f$effective_pixel_area <- f$pixel_count / f$shape_complexity
  f
}

#' @rdname compute_features
#' @export
feature_image <- function(seg, features, field, ignore_label_0 = TRUE) {

  .assert_single_layer(seg)
  .check_vector(field, "character")
  .check_columns(features, required_columns = c("id", field))
  .check_vector(ignore_label_0, "logical", 1)

  feature_mat <- as.matrix(
    features[, field, drop = FALSE]
  )

  values <- feature_image_cpp(
    labels = as.integer(terra::values(seg, mat = FALSE)),
    ids = features$id,
    feature = feature_mat,
    ignore_label_0 = ignore_label_0
  )

  out <- terra::rast(seg, nlyrs = length(field))

  terra::values(out) <- values
  names(out) <- "Feature image"
  names(out) <- field

  out
}
