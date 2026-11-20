#' Compute mean
#'
#' @description
#' Compute the mean value of raster pixels within each segment. This is a
#' lightweight alternative to [compute_features()] when only segment means are
#' required.
#'
#' @param r numeric [terra::SpatRaster-class] with one layer. Pixel values used
#'   to compute object mean.
#' @param seg single-layer [terra::SpatRaster-class]. Segmentation map
#'   of `r`, typically created with functions such as [healpix_segmentation()],
#'   [beckers_segmentation()] or [sky_grid_segmentation()], but any raster with
#'   nonnegative integer segment labels is accepted. Missing labels are excluded.
#'   A logical mask (see [direction_image()] for examples) is also accepted:
#'   `TRUE` is label 1 and `FALSE` is label 0. All selected pixels belong to
#'   one segment, even when disconnected. Use [label_components()] if connected
#'   regions should have separate labels.
#' @param return character of length one. Either `"raster"` (default) or
#'   `"vector"`, controlling whether to return a map with per-segment values or
#'   a named vector (one value per segment).
#' @param ignore_label_0 logical of length one. If `TRUE`, the segment labeled
#'   `0` is ignored.
#'
#' @details
#' Every nonmissing pixel contributes equally. Missing measurements are omitted;
#' a segment with no valid measurements receives `NA`. Included measurements
#' must otherwise be finite. For area-weighted statistics use [compute_features()].
#'
#' @return If `return = "raster"`, a [terra::SpatRaster-class] where each pixel
#'   holds its segment's feature value. If `return = "vector"`, a named numeric
#'   vector with one value per included segment, ordered by label, including
#'   segments with no valid measurements. Excluded raster cells are `NA`.
#'
#' @seealso [aggregate_by_segment()] for observations stored in a table.
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
#' compute_mean(r, seg, return = "vector")
#' plot(compute_mean(r, seg))
compute_mean <- function(r, seg,
                         return = "raster",
                         ignore_label_0 = TRUE) {

  .assert_single_layer(r)
  .assert_same_geom(r, seg)
  .assert_choice(return, c("vector", "raster"))
  labels <- .feature_labels(seg, ignore_label_0)
  values <- terra::values(r, mat = FALSE)
  .check_vector(as.numeric(values[!is.na(labels$index)]), "numeric",
                allow_na = TRUE, name = "r")
  out <- compute_mean_cpp(values, labels$index, FALSE)
  means <- out$mean[match(seq_along(labels$ids), out$id)]
  if (return == "vector") {
    return(stats::setNames(means, format(labels$ids, scientific = FALSE, trim = TRUE)))
  }
  feature_image(seg, data.frame(id = labels$ids, mean = means), "mean",
                ignore_label_0 = ignore_label_0)
}
