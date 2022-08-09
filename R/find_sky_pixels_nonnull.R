#' Find sky pixels following the non-null criteria
#'
#' Find sky pixels using the increase in the number of cells having no sky
#' pixels (the so-called null cells) as stopping criteria.
#'
#' The arguments \code{sky} and \code{slope} are passed to
#' \code{\link{thr_image}}, which output is in turn passed to
#' \code{\link{apply_thr}} along with \code{r}. As a result, \code{r} is
#' binarized and used along with \code{g} to compute the number of null cells.
#' The process is repeated but increasing \code{slope} in steps of 0.05 as long
#' as the number of null cells remains constant.
#'
#'
#' @inheritParams find_sky_pixels
#' @param slope Numeric vector of length one. See section Details in
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
#' @family Binarization Functions
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
#' bin <- find_sky_pixels_nonnull(r, sky_cs, g)
#' plot(bin)
#' }
find_sky_pixels_nonnull <- function(r, sky, g, slope = 0.5) {
  sky[is.na(sky)] <- 1
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
