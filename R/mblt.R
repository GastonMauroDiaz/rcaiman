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
#' \code{\link{fit_trend_surface}} address that topic.
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
#' which is very simple to follow.
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


#' Out-of-the-box model-based local thresholding
#'
#' Out-of-the-box version of the model-based local thresholding (MBLT)
#' algorithm.
#'
#' This function is a hard-coded version of a MBLT pipeline that starts with a
#' working binarized image and ended with a refined binarized image. The
#' pipeline combines \code{\link{find_sky_dns}},
#' \code{\link{fit_cone_shaped_model}}, \code{\link{fit_trend_surface}}, and
#' \code{\link{thr_image}}. The code can be easily inspected by calling
#' \code{ootb_mblt} --no parenthesis. Advanced users could use that code as a
#' template.
#'
#' The MBLT algorithm was first presented in
#' \insertCite{Diaz2018;textual}{rcaiman}. This version differs from that in the
#' following main aspects:
#'
#' \itemize{
#'
#' \item \eqn{intercept} is set to 0, \eqn{slope} to 1, and \eqn{w} to 0.5
#'
#' \item This version implements a regional threholding approach as first step
#' instead of a global one. Please refer to \code{\link{find_sky_dns}}. The
#' minimum number of samples (sky DNs) required is equals to the 30 % of the
#' population, considering that it is made of  \eqn{5 \times 5} sky grid cells.
#'
#' \item It does not use asynchronous acquisition under the open sky. So, the
#' cone shaped model (\code{\link{fit_cone_shaped_model}}) run without a filling
#' source, but the result of it is used as filling source for trend surface
#' fitting (\code{\link{fit_trend_surface}}).
#'
#' }
#'
#' This function search for black objects again a light background. When regular
#' canopy hemispherical image are provided as input, the algorithm will find
#' dark canopy elements against a bright sky almost everywhere in the picture.
#' If an hemispherical photograph taken under the open sky is provided, this
#' algorithm would be still searching black objects against a light background,
#' so the darker portions of the sky will be taken as objects, i.e., canopy.
#' This kind of error could be find in photographs of open forests for the same
#' reason.
#'
#' @inheritParams fit_cone_shaped_model
#'
#' @export
#' @family mblt functions
#'
#' @return Object of class list with the binarized image (named "bin") and the
#'   reconstructed skies (named "sky_cs" and "sky_s").
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
  .check_if_r_z_and_a_are_ok(r, z, a)

  bin <- find_sky_dns(r, z, a, round((360/5) * (90/5) * 0.3))

  sky_cs <- fit_cone_shaped_model(r, z, a, bin,
                                  prob = 0.95,
                                  filling_source = NULL,
                                  use_azimuth_angle = TRUE,
                                  parallel = TRUE,
                                  free_cores = 0)$image
  suppressWarnings(bin <- apply_thr(r, thr_image(sky_cs, 0, 0.5)))

  sky_s <- fit_trend_surface(r, bin,
                             m = NULL,
                             filling_source = sky_cs,
                             prob = 0.95,
                             fact = 5,
                             np = 6)$image
  suppressWarnings(bin <- apply_thr(r, thr_image(sky_s, 0, 0.5)))
  dns_8bit <- 1:15
  min_z <- Map(function(x) {
    m <- sky_s < x/255
    min(z[m], na.rm = TRUE)
  }, dns_8bit) %>% unlist()
  dns_8bit <- dns_8bit[min_z > 75]
  bin[sky_s < dns_8bit[length(dns_8bit)]/255] <- 0

  list(bin = bin, sky_cs = sky_cs, sky_s = sky_s)
}



