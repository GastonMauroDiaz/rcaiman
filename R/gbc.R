#' Gamma back correction
#'
#' Gamma back correction of JPEG images
#'
#' Digital cameras usually use sRGB as color space. It is a standard
#' developed to ensure accurate color and tone management. The transfer function
#' of sRGB, known as gamma correction, is very close to a power function with
#' the exponent 1/2.2. This is why a DN of a born-digital photograph that was
#' encoded in sRGB has a non-linear relationship with luminance despite
#' having the sensor a linear response to it.
#'
#'
#' @param DN_from_JPEG Numeric vector or object of the [SpatRaster-class]
#'   class. Digital numbers from a JPEG file (0 to 255, i.e., the standard 8-bit
#'   encoded).
#' @param gamma Numeric vector of length one. Gamma value. Please see
#'   \insertCite{Diaz2018;textual}{rcaiman} for details.
#'
#' @return The same class as the `DN_from_JPEG` argument, with dimension
#'   unchanged but values rescaled between `0` and `1` in a non-linear fashion.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(2132/2,  c(0.7836, 0.1512, -0.1558))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)/2
#'
#' caim <- expand_noncircular(caim, z, zenith_colrow)
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#'
#' bin <- apply_thr(caim$Blue, thr_isodata(caim$Blue[m]))
#'
#' display_caim(caim$Blue, bin)
#'
#' caim <- gbc(caim, 2.2)
#' }
gbc <- function(DN_from_JPEG, gamma = 2.2) {
  stopifnot(length(gamma) == 1)
  (DN_from_JPEG / 255)^gamma
}
