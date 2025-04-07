#' Normalize data using min-max rescaling
#'
#' Normalize numeric and raster data using the min-max method.
#'
#' Normalize data laying between `mn` and `mx` to the range `0`
#' to `1`. Data greater than `mx` get values greater than `1` in
#' a proportional fashion. Conversely, data less than `mn` get values less
#' than `0`.This function can be used for linear stretching of the
#' histogram.
#'
#' @param r [SpatRaster-class] or numeric vector.
#' @param mn Numeric vector of length one. Minimum expected value. Default is
#'   equivalent to enter the minimum value from `r`.
#' @param mx Numeric vector of length one. Maximum expected value. Default is
#'   equivalent to enter the maximum value from `r`.
#' @param force_range Logical vector of length one. If it is `TRUE`, the
#'   range is forced to be between `0` and `1` by flattening values
#'   found below and above those limits.
#'
#' @export
#'
#' @return An object from the same class as `r` with values from `r`
#'   linearly rescaled to make `mn` equal to zero and `mx` equal to
#'   one. Therefore, if `mn` and `mx` do not match the actual minimum
#'   and maximum from `r`, then the output will not cover the 0-to-1 range
#'   and may be outside that range if `force_range` is set to `FALSE`.
#'
#' @family Pre-processing Functions
#'
#' @examples
#' normalize_minmax(read_caim())
normalize_minmax <- function(r, mn = NULL, mx = NULL, force_range = FALSE) {
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
