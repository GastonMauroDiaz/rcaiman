#' Merge skies
#'
#' Merge skies following a weighting average approach.
#'
#' This function helps to carry on the workflow proposed in
#' \insertCite{Lang2013}{rcaiman} but automatic.
#'
#' @inheritParams ootb_mblt
#' @param sky_1,sky_2 An object of class \linkS4class{SpatRaster} produced with
#'   \code{\link{fit_coneshaped_model}}, \code{\link{fit_trend_surface}},
#'   \code{\link{fit_cie_sky_model}}, or \code{\link{ootb_sky_reconstruction}}.
#'
#' @return An object of class \linkS4class{SpatRaster} that is the result of
#'   optimally merging \code{sky_1} and \code{sky_2}.
#' @noRd
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' r <- gbc(caim$Blue)
#' r[is.na(z)] <- 0 # because FOV > 180
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- ootb_mblt(r, z, a)$bin
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' sky_points <- extract_sky_points(r, bin, g)
#' rl <- extract_rl(r, z, a, sky_points)
#' model <- fit_cie_sky_model(r, z, a, rl$sky_points,
#'                            rl$zenith_dn, sun_coord,
#'                            rmse = TRUE,
#'                            general_sky_type = "Clear")
#' sky_cie <- model$relative_luminance * model$zenith_dn
#' sky_i <- interpolate_sky_points(rl$sky_points, g) * rl$zenith_dn
#' sky_i <- cover(sky_i, sky_cie)
#' sky <- merge_skies(r, sky_cie, sky_i) # it leans toward sky_i
#' plot(sky)
#' plot(sky_i - sky)
#' plot(r/sky)
#' }
merge_skies <- function(r, sky_1, sky_2, max_w = 0.3, by = 0.01) {
  total_area <- sum(!is.na(sky_1)[], na.rm = TRUE)
  .evaluate_sky <- function(sky) {
    ratio <- r / sky
    ratio[is.infinite(ratio)] <- -0.1
    m <- ratio < 0 | ratio > 1
    area_outside_expected_values <- sum(m[], na.rm = TRUE)
    w <- area_outside_expected_values / total_area
    (sum(ratio[]^2, na.rm = TRUE) * w)
  }
  w <- seq(0, max_w, by)
  .fun <- function(i) .evaluate_sky(sky_2 * w[i] + sky_1 * (1-w[i]))
  i <- which.min(Map(.fun, seq_along(w)) %>% unlist())
  sky_2 * w[i] + sky_1 * (1-w[i])
}
