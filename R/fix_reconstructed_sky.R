#' Fix reconstructed sky
#'
#' Automatically edit a raster image of sky digital numbers (DNs) reconstructed
#' with functions such as [fit_coneshaped_model()] and
#' [fit_trend_surface()].
#'
#' The predicted sky DNs are usually erroneous near the horizon because either
#' they are a misleading extrapolation or are based on corrupted data (non-pure
#' sky DNs).
#'
#' The proposed automatic edition consists of
#' * flattening the values below
#' the minimum value from the data source defined by `r` and
#' `bin`and
#' * forcing the values toward the horizon to become gradually
#' the median value from the data source.
#'
#' The latter is achieved by calculating
#' the weighted average of the median value and the predicted sky DNs, using the
#' ratio of `z` to `90` to determine the weights.
#'
#' @param sky [SpatRaster-class]. Sky DNs predicted with functions such as
#'   [fit_coneshaped_model()] and [fit_trend_surface()].
#' @inheritParams ootb_mblt
#' @param r [SpatRaster-class]. The source of the sky DNs used to build
#'   `sky` (the data source).
#' @param bin [SpatRaster-class]. The binarization of `r` used to
#'   select the sky DNs for building the `sky` argument.
#'
#' @family Sky Reconstruction Functions
#'
#' @export
#' @return An object of class [SpatRaster-class]. The argument `sky`
#'   with dimensions unchanged but values edited.
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#'
#'
#' sky_points <- extract_rl(r, z, a, extract_sky_points_simple(r, z, a),
#'                          NULL,
#'                          use_window = FALSE)# this is important when
#'                                             # extract_sky_points_simple()
#'                                             # is used
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' summary(model$model)
#' sky_cs <- model$fun(z, a)
#' plot(r/sky_cs)
#' plot(sky_cs)
#' persp(terra::aggregate(sky_cs, 20), theta = 90, phi = 20)
#' bin <- apply_thr(r, thr_mblt(sky_cs, 0, 0.5))
#' plot(bin)
#' sky_cs <- fix_reconstructed_sky(sky_cs, z, r, bin)
#' bin <- apply_thr(r, thr_mblt(sky_cs, 0, 0.5))
#' plot(bin)
#' persp(terra::aggregate(sky_cs, 20), theta = 90, phi = 20)
#' }
fix_reconstructed_sky <- function(sky, z, r, bin) {
  .is_single_layer_raster(r, "r")
  .was_normalized(r)
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(sky, "sky")
  stopifnot(.get_max(z) <= 90)
  .is_single_layer_raster(bin, "bin")
  terra::compareGeom(r, z)
  terra::compareGeom(r, sky)
  terra::compareGeom(r, bin)
  .is_logic_and_NA_free(bin, "bin")

  mn <- min(r[bin])
  sky[sky < mn] <- mn
  x <- median(r[bin])
  w <- z / 90
  x * w + sky * (1 - w)
}



