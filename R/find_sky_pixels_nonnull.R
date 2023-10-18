#' Find sky pixels following the non-null criteria
#'
#' Cells without sky pixels are the so-called null cells. This type of cells are
#' mathematically intractable by models typically used to obtain canopy metrics.
#' This function find sky pixels using increase in number of null cells as the
#' stopping criteria.
#'
#' The arguments `sky`, `intercept`, `slope`, and `w` are passed to [thr_mblt()]
#' whose output is in turn passed to [apply_thr()] along with `r`. As a result,
#' `r` is binarized and used along with `g` to compute the number of null cells.
#' The process is repeated but increasing `w` in steps of 0.05 as long as
#' the number of null cells remains constant.
#'
#'
#' @inheritParams find_sky_pixels
#' @inheritParams thr_mblt
#' @inheritParams ootb_mblt
#' @inheritParams extract_sky_points
#' @param sky An object of class [SpatRaster-class] produced with
#'   [fit_coneshaped_model()], [fit_trend_surface()], [fit_cie_sky_model()], or
#'   [ootb_sky_reconstruction()]. It also support a numeric vector of length
#'   one. For instance, it could be a value obtained with a combination of
#'   [extract_sky_points()] and [extract_dn()]. The latter can be understood as
#'   modelling the sky with a plane.
#'
#' @return An object of class [SpatRaster-class] with values `0` and `1`.
#' @export
#'
#' @family Binarization Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize(., 0, 2^16)
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' bin <- ootb_obia(caim, z, a, gamma = NULL)
#' g <- sky_grid_segmentation(z, a, 3)
#' r <- caim$Blue
#' sky_points <- extract_sky_points(r, bin, g,
#'                                  dist_to_plant = 5,
#'                                  min_raster_dist = 5)
#' rl <- extract_rl(r, z, a, sky_points)
#' model <- fit_coneshaped_model(rl$sky_points)
#' summary(model$model)
#'
#' sky <- model$fun(z, a)
#' sky <- fit_trend_surface(sky, z, a, !is.na(z))$image
#' plot(r/sky)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels_nonnull(r, sky, g)
#' plot(bin)
#' }
find_sky_pixels_nonnull <- function(r, sky, g,
                                    intercept = 0,
                                    slope = 1,
                                    w = 0.5) {
  sky[is.na(sky)] <- 1
  .get_nulls_no <- function(w) {
    thr <- suppressWarnings(thr_mblt(sky, intercept, slope* w))
    bin <- apply_thr(r, thr)
    no_of_nulls <- extract_feature(bin, g, mean, return_raster = FALSE)
    sum(no_of_nulls == 0)
  }
  unlock <- TRUE
  while (unlock) {
    no_nulls_1 <- .get_nulls_no(w)
    w <- w + 0.05
    no_nulls_2 <- .get_nulls_no(w)
    unlock <- no_nulls_2 <= no_nulls_1
    if (!unlock) w <- w - 0.05
  }
  thr <- suppressWarnings(thr_mblt(sky, intercept, slope * w))
  apply_thr(r, thr)
}
