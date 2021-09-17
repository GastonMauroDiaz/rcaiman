
#' Normalize data
#'
#' Normalize data lying between \code{mn} and \code{mx} in the range \code{0} to
#' \code{1}. Data greater than \code{mx} get values greater than \code{1} in a
#' proportional fashion. Conversely, data less than \code{mn} get values less
#' than \code{0}.
#'
#' @param r \code{\linkS4class{Raster}} or numeric vector.
#' @param mn Numeric vector of length one. Minimum expected value.
#' @param mx Numeric vector of length one. Maximum expected value.
#'
#' @examples
#' normalize(read_caim(), 0, 255)
normalize <- function(r, mn, mx) {
  stopifnot(length(mn) == 1)
  stopifnot(length(mx) == 1)
  (r - mn) / (mx - mn)
}

#' Gamma back correction
#'
#' @param DN_from_JPEG Numeric vector or an object from the
#'   \code{\linkS4class{Raster}} class. Digital numbers from a JPEG file.
#' @param gamma Numeric vector of length one. Gamma value. Please see the
#'   reference for details.
#'
#' @return Normalized \code{\linkS4class{Raster}}.
#' @export
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
