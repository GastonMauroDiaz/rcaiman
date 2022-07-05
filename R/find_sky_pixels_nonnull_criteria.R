#' Find sky pixels following the non-null criteria
#'
#' To produce a binarized image, the arguments \code{sky} and \code{slope} is
#' passed to \code{\link{thr_image}}, which result is in turn passed to
#' \code{\link{apply_thr}} along with \code{r}.
#'
#' A sky grid of \eqn{10 \times 10} degrees is used to compute the number of
#' cells having none sky pixels (the so-called null cells). The process is
#' repeated but increasing \code{slope} in steps of 0.05 as long as the number
#' of null cells remain constant.
#'
#'
#' @inheritParams find_sky_pixels
#' @param slope Numeric vector of length one. Please, see the Details section of
#'   \code{\link{thr_image}}.
#' @inheritParams extract_sky_points
#' @param sky An object of class \linkS4class{SpatRaster} produced with
#'   \code{\link{fit_coneshaped_model}}, \code{\link{fit_trend_surface}},
#'   \code{\link{fit_cie_sky_model}}, or \code{\link{ootb_sky_reconstruction}}.
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' bin <- find_sky_pixels(r, z, a)
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_rl(r, z, a, sky_points, NULL)
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' sky_cs <- model$fun(z, a)
#' g[mask_hs(z, 0, 10) | mask_hs(z, 70, 90)] <- NA
#' bin <- find_sky_pixels_nonnull_criteria(r, sky_cs, g)
#' plot(bin)
#' }
find_sky_pixels_nonnull_criteria <- function(r, sky, g, slope = 0.5) {
  .get_nulls_no <- function(slope) {
    thr <- suppressWarnings(thr_image(sky, 0, slope))
    bin <- apply_thr(r, thr)
    no_of_nulls <- extract_feature(bin, g, mean, return_raster = FALSE)
    sum(no_of_nulls == 0)
  }
  unlock <- TRUE
  while (unlock) {
    no_nulls_1 <- .get_nulls_no(slope)
    slope <- slope + 0.05
    no_nulls_2 <- .get_nulls_no(slope)
    unlock <- no_nulls_2 <= no_nulls_1
    if (!unlock) slope <- slope - 0.05
  }
  thr <- suppressWarnings(thr_image(sky, 0, slope))
  apply_thr(r, thr)
}
