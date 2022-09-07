#' Out-of-the-box object-based image analysis of canopy photographs
#'
#' Out-of-the-box version of methods first presented in
#' \insertCite{Diaz2015;textual}{rcaiman}.
#'
#' This function is a hard-coded version of a pipeline that combines these main
#' functions \code{\link{enhance_caim}}, \code{\link{polar_qtree}}, and
#' \code{\link{obia}}. The code can be easily inspected by calling
#' \code{ootb_obia} --no parenthesis. Advanced users can use that code as a
#' template.
#'
#' Pixels from the synthetic layer returned by \code{\link{obia}} that lay
#' between \code{0} and \code{1} are assigned to the class \emph{plant} only if
#' they are:
#'
#' \itemize{
#'
#' \item \code{0} after
#' \code{\link{defuzzify}} with a sky grid segmentation of \code{10} degrees,
#'
#' \item \code{0} after \code{\link{apply_thr}} with a threshold of
#' \code{0.5}, and
#'
#' \item not exclusively surrounded by sky pixels.
#'
#' }
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
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(2132/2, lens("Olloclip"))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)/2
#' caim <- expand_noncircular(caim, z, zenith_colrow) %>% normalize()
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#' bin <- ootb_obia(caim, z, a)
#' plot(bin)
#' }
ootb_obia <- function(caim, z, a, m = NULL, sky_blue = NULL) {
  ecaim <- enhance_caim(caim, m, sky_blue = sky_blue,
                        w_red = 0, gamma = 2.2, thr = 0.5,
                        fuzziness = 100)
  if (is.null(m)) m <- !is.na(z)
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
  r <- gbc(caim$Blue*255, gamma = 2.2)
  synth <- obia(r, z, a, bin, seg)
  synth <- terra::cover(synth, bin)
  synth[!bin] <- 0

  g <- sky_grid_segmentation(z, a, 10)
  bin_obia <- defuzzify(synth, g) | apply_thr(synth, 0.5)
  ma <- matrix(c(1,1,1,1,-8,1,1,1,1), ncol = 3, nrow = 3)
  bin_obia[terra::focal(bin_obia, ma) == 8] <- 1
  bin_obia[!bin] <- 0
  as.logical(bin_obia)
}
