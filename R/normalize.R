#' Normalize data
#'
#' Normalize numeric and raster data.
#'
#' Normalize data laying between \code{mn} and \code{mx} to the range \code{0}
#' to \code{1}. Data greater than \code{mx} get values greater than \code{1} in
#' a proportional fashion. Conversely, data less than \code{mn} get values less
#' than \code{0}.This function can be used for linear stretching of the
#' histogram.
#'
#' @param r \code{\linkS4class{SpatRaster}} or numeric vector.
#' @param mn Numeric vector of length one. Minimum expected value. Default is
#'   equivalent to enter the minimum value from \code{r}.
#' @param mx Numeric vector of length one. Maximum expected value. Default is
#'   equivalent to enter the maximum value from \code{r}.
#' @param force_range Logical vector of length one. If it is \code{TRUE}, the
#'   range is forced to be between \code{0} and \code{1} by flattening values
#'   found below and above those limits.
#'
#' @export
#'
#' @return An object from the same class as \code{r} with values from \code{r}
#'   linearly rescaled to make \code{mn} equal to zero and \code{mx} equal to
#'   one. Therefore, if \code{mn} and \code{mx} do not match the actual minimum
#'   and maximum from \code{r}, then the output will not cover the 0-to-1 range
#'   and may be outside that range if \code{force_range} is set to \code{FALSE}.
#'
#' @family Pre-processing Functions
#'
#' @examples
#' normalize(read_caim(), 0, 255)
normalize <- function(r, mn = NULL, mx = NULL, force_range = FALSE) {
  if (is.null(mn)) mn <- .get_min(r)
  if (is.null(mx)) mx <- .get_max(r)
  stopifnot(length(mn) == 1)
  stopifnot(length(mx) == 1)
  r <- (r - mn) / (mx - mn)
  if (force_range) {
    r[r < 0] <- 0
    r[r > 1] <- 1
  }
  r
}
