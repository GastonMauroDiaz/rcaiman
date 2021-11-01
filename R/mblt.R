#' Threshold image
#'
#' Transform background digital number into threshold values.
#'
#' Transform background digital number into threshold values by means of the
#' Equation 1 presented in \insertCite{Diaz2018;textual}{rcaiman}, which is a
#' linear function with the slope modified by a weighting parameter. This simple
#' function was found by studying canopy models, also known as targets, which
#' are planes with holes made of a rigid and dark material. These models were
#' backlighted with homogeneous lighting, photographed with a Nikon Coolpix 5700
#' set to acquire in JPEG format, and those images were gamma back corrected
#' with a default gamma value equal to 2.2 (see \code{\link{gbc}}). Results
#' clearly shown that the optimal threshold value was linearly related with the
#' background digital number. Therefore, that shifts the aim from finding the
#' optimal threshold to obtaining the background DN as if the canopy were not
#' there. Functions \code{\link{fit_cone_shaped_model}} and
#' \code{\link{fit_trend_surface_to_sky_dn}} address that topic.
#'
#' It is worth noting that Equation 1 was developed with 8-bit images, so
#' calibration of new coefficient should be done in the 0 to 255 domain since
#' that is what \code{\link{thr_image}} expect, although the input \code{dn}
#' should be normalized. The latter --that might sound counter intuitive-- was a
#' design decision aiming to harmonize the whole package.
#'
#' To apply the weighting parameter (w) from Equation 1, just provide the
#' argument \code{slope} as \code{slope_value * w}.
#'
#' Type \code{thr_image} --no parenthesis-- in the console to inspect the code,
#' which is very simple.
#'
#' @param dn Numeric vector or \linkS4class{RasterLayer}. Digital number of the
#'   background. These values should be normalized and, if they are extracted
#'   from JPEG image, gamma back corrected.
#' @param intercept,slope Numeric vector of length one. These are linear
#'   function coefficients. Please, see the Details section of
#'   \code{\link{thr_image}}.
#'
#' @export
#'
#' @family mblt functions
#' @seealso \code{\link{normalize}}, \code{\link{gbc}}, \code{\link{apply_thr}}
#'   and \code{\link{regional_thresholding}}.
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' thr_image(gbc(125), -8, 1)
thr_image <- function (dn, intercept, slope) {
  stopifnot(length(intercept) == 1)
  stopifnot(class(intercept) == "numeric")
  stopifnot(length(slope) == 1)
  stopifnot(class(slope) == "numeric")

  dn <- dn * 255

  if (.get_max(dn) > 255) warning("\"dn\" values should be normalized")
  dn[dn > 255] <- 255

  thr <- intercept  + slope * dn
  thr[thr < 0] <- 0

  thr / 255
}


#' Out-of-the-box Model-based local thresholding
#'
#' Out-of-the-box version of the model-based local thresholding (MBLT)
#' algorithm.
#'
#' This function is a hard-coded version of a MBLT pipeline that starts with a
#' working binarized image and ended with a refined binarized image. The
#' pipeline combines \code{\link{find_sky_dns}},
#' \code{\link{fit_cone_shaped_model}},
#' \code{\link{fit_trend_surface_to_sky_dn}}, and \code{\link{thr_image}}. The
#' code can be easily inspected by calling \code{ootb_mblt} --no parenthesis.
#' Advanced users could use that code as a template.
#'
#' The MBLT algorithm was first presented in
#' \insertCite{Diaz2018;textual}{rcaiman}. This version differs from that in the
#' following main aspects:
#'
#' \itemize{
#'
#' \item $intercept$ is set to 0, $slope$ to 1, and $w$ to 0.5
#'
#' \item This version implements a regional threholding approach as first step
#' instead of a global one. Please refer to \code{\link{find_sky_dns}}. The
#' minimum number of samples (sky dns) required is equals to the 30 % of the
#' population, considering that it is made of 5 $/times$ 5 sky grid cells.
#'
#' \item It does not use asynchronous acquisition under the open sky. So, the
#' cone shaped model (\code{\link{fit_cone_shaped_model}}) run without a
#' filling source, but the result of it is used as filling source for trend
#' surface fitting (\code{\link{fit_trend_surface_to_sky_dn}}).
#'
#' \item The sDN obtained by trend surface fitting is merged with the sDN
#' obtained with \code{\link{fit_cone_shaped_model}}. To merge them, a weighted
#' average is calculated, being weights calculated as $/theta^2 / 90^2$ (Near
#' the zenith, values obtained by means of trend surface fitting prevail over
#' the ones obtained with the cone shaped model, and the opposite occur near the
#' horizon).
#'
#' }
#'
#'
#' @inheritParams fit_cone_shaped_model
#'
#' @export
#' @family mblt functions
#'
#' @return \link{RasterStack} with the binarized image (layer named "bin") and
#'   the reconstructed skies. The layer named "sky_m" is the cone shaped model,
#'   "sky_s" is the trend surface, and "sky" is the combination of both, and
#'   actually used to obtain "bin".
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' bin <- ootb_mblt(blue, z, a)$bin
#' plot(bin)
#' }
ootb_mblt <- function(r, z, a) {
  bin <- find_sky_dns(r, z, a, round((360/5) * (90/5) * 0.3))

  sky_m <- fit_cone_shaped_model(r, z, a, bin)$image
  bin <- apply_thr(r, thr_image(sky_m, 0, 0.5))

  m <- mask_image(z)
  sky_s <- fit_trend_surface_to_sky_dn(r, m, bin, sky_m)$image
  w <- z^2 / 90^2
  sky <- sky_s * (1 - w) + sky_m *  w
  bin <- apply_thr(r, thr_image(sky, 0, 0.5))

  r <- stack(bin, sky_m, sky_s, sky)
  names(r) <- c("bin", "sky_m", "sky_s", "sky")
  r
}






