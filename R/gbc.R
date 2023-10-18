#' Gamma back correction
#'
#' Gamma back correction of JPEG images
#'
#' Digital cameras usually use sRGB as color space, which is a standard
#' developed to ensure accurate color and tone management. The transfer function
#' of sRGB, known as gamma correction, is very close to a power function with
#' the exponent 1/2.2. This is why a DN of a born-digital photograph that was
#' encoded in sRGB has a non-linear relationship with luminance despite
#' having the sensor a linear response.
#'
#'
#' @param DN_from_JPEG Numeric vector or object from the [SpatRaster-class]
#'   class. Digital numbers from a JPEG file (0 to 255, i.e., the standard 8-bit
#'   encoded).
#' @param gamma Numeric vector of length one. Gamma value. Please see
#'   \insertCite{Diaz2018;textual}{rcaiman} for details.
#'
#' @return The same class as `DN_from_JPEG`, with dimension unchanged but values
#'   rescaled between `0` and `1` in a non-linear fashion.
#' @export
#'
#' @family Pre-processing Functions
#'
#'
#' @references \insertAllCited{}
#'
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' r
#' gbc(r)
gbc <- function(DN_from_JPEG, gamma = 2.2) {
  stopifnot(length(gamma) == 1)
  (DN_from_JPEG / 255)^gamma
}
