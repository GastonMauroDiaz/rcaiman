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
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path), z, a)
#'
#' sky <- ootb_build_sky_vault(r, z, a, ootb_sky$sky_points, ootb_sky)
#'
#' calc_oor_index(r, sky$sky)
#' }
calc_oor_index <- function(r, sky) {
  ratio <- r / sky
  ratio[is.infinite(ratio)] <- 1e+10
  out.of.range_ratio <- ratio - normalize_minmax(ratio, 0, 1, TRUE)
  out.of.range_ratio <- sum(out.of.range_ratio[]^2,
                            na.rm = TRUE)
  out.of.range_ratio
}
