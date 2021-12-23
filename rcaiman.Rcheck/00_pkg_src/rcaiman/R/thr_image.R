#' Threshold image
#'
#' Transform background digital number into threshold values.
#'
#' This function transforms background digital number into threshold values by
#' means of the Equation 1 presented in \insertCite{Diaz2018;textual}{rcaiman},
#' which is a linear function with the slope modified by a weighting parameter.
#' This simple function was found by studying canopy models, also known as
#' targets, which are planes with holes made of a rigid and dark material. These
#' models were backlighted with homogeneous lighting, photographed with a Nikon
#' Coolpix 5700 set to acquire in JPEG format, and those images were gamma back
#' corrected with a default gamma value equal to 2.2 (see \code{\link{gbc}}).
#' Results clearly shown that the optimal threshold value was linearly related
#' with the background digital number. Therefore, that shifts the aim from
#' finding the optimal threshold to obtaining the background DN as if the canopy
#' was not there. Functions \code{\link{fit_cone_shaped_model}} and
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
#' @family MBLT functions
#' @seealso \code{\link{normalize}}, \code{\link{gbc}}, \code{\link{apply_thr}}
#'   and \code{\link{regional_thresholding}}.
#'
#' @references \insertAllCited{}
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
