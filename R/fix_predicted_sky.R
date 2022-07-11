#' Fix predicted sky
#'
#' Automatically edit a raster image of sky digital numbers (DNs) predicted with
#' functions such as \code{\link{fit_coneshaped_model}} and
#' \code{\link{fit_trend_surface}}.
#'
#' The predicted sky DNs are usually erroneous near the horizon because they are
#' a misleading extrapolation or are based on corrupted data (non-pure sky DNs).
#'
#' This automatic edition consists of (1) flattening values below the minimum
#' input data (2) and forcing the values toward the horizon to gradually become
#' the median input data. The latter is achieved by calculating the weighted
#' average of the median value and the predicted sky DNs --\code{z / 90} is used
#' to determine the weights.
#'
#' @param sky \linkS4class{SpatRaster}. Sky DNs predicted with functions such
#'   as \code{\link{fit_coneshaped_model}} and \code{\link{fit_trend_surface}}.
#' @param z \linkS4class{SpatRaster}. The result of a call to
#'   \code{\link{zenith_image}}.
#' @param r \linkS4class{SpatRaster}. The source of the sky DNs used to build
#'   \code{sky}.
#' @param bin \linkS4class{SpatRaster}. The binarization of \code{r} used to
#'   select the sky DNs for building \code{sky}.
#'
#' @family Sky reconstruction functions
#'
#' @export
#' @return An object of class \linkS4class{SpatRaster}. The edited
#'   version of \code{sky}. Pixel dimensions should remain unchanged.
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels(r, z, a)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_rl(r, z, a, sky_points, NULL)
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' sky_cs <- model$fun(z, a)
#' sky_cs <- fix_predicted_sky(sky_cs, z, r, bin)
#' persp(sky_cs, theta = 90, phi = 0)
#' }
fix_predicted_sky <- function(sky, z, r, bin) {
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
  x <- quantile(r[bin], 0.5)
  w <- z / 90
  x * w + sky * (1 - w)
}
