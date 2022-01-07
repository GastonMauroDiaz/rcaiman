#' Out-of-the-box model-based local thresholding
#'
#' Out-of-the-box version of the model-based local thresholding (MBLT)
#' algorithm.
#'
#' This function is a hard-coded version of a MBLT pipeline that starts with a
#' working binarized image and ends with a refined binarized image. The pipeline
#' combines \code{\link{find_sky_dns}}, \code{\link{fit_cone_shaped_model}},
#' \code{\link{fit_trend_surface}}, and \code{\link{thr_image}}. The code can be
#' easily inspected by calling \code{ootb_mblt} --no parenthesis. Advanced users
#' can use that code as a template.
#'
#' The MBLT algorithm was first presented in
#' \insertCite{Diaz2018;textual}{rcaiman}. The version presented here differs
#' from that in the following main aspects:
#'
#' \itemize{
#'
#' \item \eqn{intercept} is set to 0, \eqn{slope} to 1, and \eqn{w} to 0.5
#'
#' \item This version implements a regional threholding approach as first step
#' instead of a global one. Please refer to \code{\link{find_sky_dns}}. The
#' minimum number of samples (sky DNs) required is equals to the 30 percent of
#' the population, considering that it is made of  \eqn{5 \times 5} sky grid
#' cells.
#'
#' \item It does not use asynchronous acquisition under the open sky. So, the
#' cone shaped model (\code{\link{fit_cone_shaped_model}}) run without a filling
#' source, but the result of it is used as filling source for trend surface
#' fitting (\code{\link{fit_trend_surface}}).
#'
#' \item If the cone shaped model predicts values below zero, those values are
#' set to zero and values toward the horizon are forced to gradually become the
#' median sky DN calculated from the DNs finded by \code{\link{find_sky_dns}}.
#'
#' }
#'
#' This function searches for black objects against a light background. When
#' regular canopy hemispherical images are provided as input, the algorithm will
#' find dark canopy elements against a bright sky almost everywhere in the
#' picture, and the result will fit user expectations. However, if an
#' hemispherical photograph taken under the open sky is provided, this algorithm
#' would be still searching black objects against a light background, so the
#' darker portions of the sky will be taken as objects, i.e., canopy. As a
#' consequence, this will not fit users expectations, since they require the
#' classes ‘Gap’ and ‘No-gap’. This kind of error could be find in photographs
#' of open forests for the same reason.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018}{rcaiman}.
#'
#' @inheritParams fit_cone_shaped_model
#'
#' @export
#' @family MBLT functions
#'
#' @return Object of class list with the binarized image (named ‘bin’) and the
#'   reconstructed skies (named ‘sky_cs’ and ‘sky_s’).
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/4_D_2_DSCN4502.JPG", package = "rcaiman")
#' r <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' bin <- ootb_mblt(blue, z, a)
#' plot(bin$bin)
#' }
ootb_mblt <- function(r, z, a) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  bin <- find_sky_dns(r, z, a, round((360/5) * (90/5) * 0.3))
  sky_cs <- fit_cone_shaped_model(r, z, a, bin,
                                  prob = 0.95,
                                  filling_source = NULL,
                                  use_azimuth_angle = TRUE,
                                  parallel = TRUE,
                                  free_cores = 0)$image
  if (min(sky_cs[], na.rm = TRUE) < 0) {
    sky_cs[sky_cs < 0] <- 0
    x <- quantile(r[bin], 0.5)
    w <- z / 90
    sky_cs <- x * w + sky_cs * (1 - w)
  }
  suppressWarnings(bin <- apply_thr(r, thr_image(sky_cs, 0, 0.5)))
  sky_s <- fit_trend_surface(r, bin,
                             m = NULL,
                             filling_source = sky_cs,
                             prob = 0.95,
                             fact = 5,
                             np = 6)$image
  suppressWarnings(bin <- apply_thr(r, thr_image(sky_s, 0, 0.5)))

  bin[is.na(z)] <- 0
  list(bin = bin, sky_cs = sky_cs, sky_s = sky_s)
}
