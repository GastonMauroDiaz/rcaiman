#' Invert gamma correction of JPEG images
#'
#' @description
#' Approximates the inversion of the gamma encoding applied to JPEG
#' images.
#'
#' @details
#' Digital cameras typically encode images using the sRGB color space, which
#' applies a non-linear transformation—commonly referred to as gamma
#' correction—to the sensor's linear luminance response
#' \insertCite{Diaz2018}{rcaiman}. This function applies a power transformation
#' to approximate the inverse of that encoding, restoring a response closer to
#' linear but with uncertainties.
#'
#' In \insertCite{Macfarlane2014;textual}{rcaiman}, page 252, we can read: “our
#' comparison of gamma-adjusted raw data and incamera JPG images suggests that
#' in-camera processing results in more substantial alteration of digital
#' numbers and that this alteration would be difficult to back-correct to truly
#' linear values.”
#'
#' @param dn numeric vector or [terra::SpatRaster-class]. Digital
#'   numbers from a JPEG file (range 0–255, as per standard 8-bit encoding).
#' @param gamma numeric vector of length one. Exponent applied in the inverse
#'   gamma correction (typically 2.2 for sRGB).
#'
#' @returns Same properties as `dn`, with values adjusted by
#'   inverse gamma correction and rescaled to the range \eqn{[0, 1]}.
#'
#' @export
#'
#' @seealso [thr_mblt()]
#'
#' @references \insertAllCited{}
#'
#' @examples
#'  invert_gamma_correction(125)
invert_gamma_correction <- function(dn, gamma = 2.2) {
  .check_vector(gamma, "numeric", 1, sign = "positive")
  handling_dn <- c(
    tryCatch(.check_vector(dn, "numeric", sign = "any"),
             error = function(e) FALSE),
    tryCatch(.assert_spatraster(dn),
             error = function(e) FALSE)
  )
  if (!any(handling_dn)) {
    stop("`dn` must be either numeric vector or SpatRaster class.")
  }
  (dn / 255)^gamma
}
