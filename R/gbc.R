#' Gamma back correction
#'
#' Gamma back correction of JPEG images.
#'
#' @param DN_from_JPEG Numeric vector or an object from the
#'   \code{\linkS4class{Raster}} class. Digital numbers from a JPEG file.
#' @param gamma Numeric vector of length one. Gamma value. Please see the
#'   reference for details.
#'
#' @return Normalized \code{\linkS4class{Raster}}.
#' @export
#'
#' @seealso \code{\link{normalize}}
#' @family Tools functions
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' r <- read_caim()
#' gbc(r)
gbc <- function(DN_from_JPEG, gamma = 2.2) {
  stopifnot(length(gamma) == 1)
  (DN_from_JPEG / 255)^gamma
}
