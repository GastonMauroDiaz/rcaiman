#' Out-of-the-box object-based image analysis of hemispherical photographs
#'
#' Out-of-the-box version of methods first presented in
#' \insertCite{Diaz2015;textual}{rcaiman}.
#'
#' This function is a hard-coded version of a pipeline that combines these main
#' functions \code{\link{enhance_caim}}, \code{\link{polar_qtree}},
#' \code{\link{obia}}. The code can be easily inspected by calling
#' \code{ootb_obia} --no parenthesis. Advanced users can use that code as a
#' template.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package.
#'
#' @inheritParams enhance_caim
#' @inheritParams ootb_mblt
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#'
#'
#' @family Binarization Functions
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' sky_blue_sample <- crop(caim, ext(610,643,760,806))
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' bin <- ootb_obia(caim, z, a, !is.na(z), sky_blue)
#' plot(bin)
#' }
ootb_obia <- function(caim, z, a, m, sky_blue) {
  ecaim <- enhance_caim(caim, m, sky_blue = sky_blue,
                        w_red = 0, gamma = 2.2, thr = 0.5,
                        fuzziness = 100)
  sunlit_canopy <- mask_sunlit_canopy(caim, m)
  ecaim[sunlit_canopy] <- 0
  mx <- quantile(ecaim[], 0.9)
  ecaim <- normalize(ecaim, 0, mx, force_range = TRUE)

  bin <- find_sky_pixels(ecaim, z, a)
  g <- sky_grid_segmentation(z, a, 10)
  sky_points <- extract_sky_points(ecaim, bin, g)
  sky_points <- extract_rl(ecaim, z, a, sky_points, NULL)
  sky_flat <- z
  sky_flat[] <- mean(sky_points$sky_points$dn)
  g[mask_hs(z, 0, 10) | mask_hs(z, 70, 90)] <- NA
  bin <- find_sky_pixels_nonnull(ecaim, sky_flat, g)

  seg <- polar_qtree(caim, z, a, scale_parameter = 0.2)
  r <- gbc(caim$Blue*255)
  obia(r, z, a, bin, seg)
}