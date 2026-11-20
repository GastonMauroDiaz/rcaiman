#' Smooth a histogram with a Gaussian kernel
#'
#' @description
#' Smooth histogram frequencies while preserving the bins, their units, and
#' the total frequency or weight.
#'
#' @details
#' A discrete Gaussian kernel extends to `ceiling(3 * sigma)` bins on each
#' side. Frequencies are reflected at the ends, without repeating the endpoint
#' bin. The smoothed frequencies are rescaled to preserve their original sum;
#' `density` is then recalculated. Reflection and smoothing can affect peaks
#' near the ends of the histogram.
#'
#' The input is not modified. Keeping both versions allows peak detection on
#' the original histogram and corner construction on the smoothed histogram,
#' as in [thr_two_corner()].
#'
#' @param h histogram from [build_histogram()], or a list with equally spaced
#'   `breaks`, matching `mids`, and nonnegative finite `counts` with positive sum.
#' @param sigma positive numeric vector of length one. Standard deviation of
#'   the Gaussian kernel, measured in bins. Default 2. `NULL` returns `h`
#'   unchanged.
#'
#' @return A histogram with the same bins and updated `counts` and `density`.
#' @seealso [find_histogram_peaks()], [thr_rosin()]
#' @export
#' @examples
#' x <- qbeta(ppoints(1000), 2, 10)
#' h <- build_histogram(x)
#' h_smooth <- smooth_histogram(h, sigma = 2)
#' plot(h_smooth, freq = FALSE)
#' lines(h$mids, h$density, col = "grey")
smooth_histogram <- function(h, sigma = 2) {
  .check_threshold_histogram(h)
  .check_vector(sigma, "numeric", 1, sign = "positive", allow_null = TRUE)
  if (is.null(sigma)) return(h)
  # Normalize before convolution, then restore the original mass.
  total <- sum(h$counts)
  counts <- .smooth_hist_gaussian(h$counts / total, sigma)
  h$counts <- counts / sum(counts) * total
  h$density <- h$counts / sum(h$counts) / diff(h$breaks)
  h
}
