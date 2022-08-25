#' Gamma back correction
#'
#' Gamma back correction of JPEG images.
#'
#' @param DN_from_JPEG Numeric vector or object from the
#'   \code{\linkS4class{SpatRaster}} class. Digital numbers from a JPEG file (0
#'   to 255, i.e., the standard 8-bit encoded).
#' @param gamma Numeric vector of length one. Gamma value. Please see
#'   \insertCite{Diaz2018;textual}{rcaiman} for details.
#'
#' @return The same class as \code{DN_from_JPEG}, with dimension unchanged but
#'   values rescaled between \code{0} and \code{1} in a non-linear fashion.
#' @export
#'
#' @family Pre-processing Functions
#'
#'
#' @references \insertAllCited{}
#'
#' @examples
#' r <- read_caim()
#' r
#' gbc(r)
gbc <- function(DN_from_JPEG, gamma = 2.2) {
  stopifnot(length(gamma) == 1)
  (DN_from_JPEG / 255)^gamma
}
