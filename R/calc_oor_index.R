#' Calc out-of-range index
#'
#' The __out-out-of-range__ (OOR) index is calculated as foollow:
#'
#' \eqn{\sum_{i = 1}^{N}(r_i/sky_i)^2},
#'
#' where \eqn{r} is a canopy image with radiance data, \eqn{sky} is an
#' image with unobstructed sky radiance data, \eqn{i} is the index that
#' represents the position of a given pixel on the raster grid, and \eqn{N} is
#' the total number of pixels that satisfy: \eqn{r_i/sky_i<0} or
#' \eqn{r_i/sky_i>1}.
#'
#' @param r [SpatRaster-class]. The blue channel of a canopy photograph.
#' @param sky [SpatRaster-class]. The blue radiance of the unobscured sky
#'
#' @returns Numeric vector of length one.
#'
#' @family Tool Functions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' r <- correct_vignetting(r, z, c(0.0638, -0.101)) %>% normalize()
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 30), "thr_isodata")
#' bin <- bin & select_sky_vault_region(z, 0, 80)
#' sky_points <- extract_sky_points(r, bin, sky_grid_segmentation(z, a, 3))
#' sky_points <- extract_rel_radiance(r, z, a, sky_points, no_of_points = NULL)
#'
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' summary(model$model)
#' sky_cs <- model$fun(z, a)
#' plot(r/sky_cs)
#' calc_oor_index(r, sky_cs)
#' }
calc_oor_index <- function(r, sky) {
  ratio <- r / sky
  ratio[is.infinite(ratio)] <- 1e+10
  out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
  out.of.range_ratio <- sum(out.of.range_ratio[]^2,
                            na.rm = TRUE)
  out.of.range_ratio
}
