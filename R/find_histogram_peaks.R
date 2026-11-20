#' Find local peaks in a histogram
#'
#' @description
#' Return the bin indices of positive local maxima in a supplied histogram.
#' Smoothing is a separate step, controlled by [smooth_histogram()].
#'
#' @details
#' A peak exceeds its immediate neighbours. A flat-topped peak contributes
#' one index at the middle of its plateau; for an even number of bins the
#' lower index is used. Endpoint peaks are included when they exceed their
#' single neighbour. A constant histogram has no local peaks. No prominence
#' filter, mode-count assumption, or significance test is applied.
#'
#' Peaks can be explored across several smoothing widths. Persistence across
#' widths describes the supplied distribution; it does not establish a
#' population's number of modes. Weighted histograms are treated in exactly
#' the same way as unweighted histograms.
#'
#' @inheritParams smooth_histogram
#' @return Integer vector of bin indices in increasing order, or `integer(0)`
#'   if no peaks are found. Use `h$mids[peaks]` for positions in data units.
#' @seealso [find_macfarlane_peaks()], [thr_rosin()], [thr_two_corner()]
#' @export
#' @examples
#' x <- c(qnorm(ppoints(1500), -1, .2), qnorm(ppoints(1000), 1, .3))
#' h <- build_histogram(x)
#' h_smooth <- smooth_histogram(h)
#' peaks <- find_histogram_peaks(h_smooth)
#' h$mids[peaks]
#' plot(h_smooth, freq = FALSE)
#' points(h_smooth$mids[peaks], h_smooth$density[peaks], pch = 19)
#'
#' # Explore how the number and positions of peaks change with smoothing.
#' sigmas <- c(.5, 1, 2, 4, 8, 16)
#' peaks_by_sigma <- lapply(sigmas, function(sigma) {
#'   find_histogram_peaks(smooth_histogram(h, sigma))
#' })
#' data.frame(sigma = sigmas, n_peaks = lengths(peaks_by_sigma))
#' positions <- lapply(peaks_by_sigma, function(i) h$mids[i])
#' plot(rep(sigmas, lengths(positions)), unlist(positions),
#'      xlab = "Smoothing sigma (bins)", ylab = "Peak position", pch = 19)
find_histogram_peaks <- function(h) {
  .check_threshold_histogram(h)
  runs <- rle(h$counts)
  if (length(runs$values) == 1L) return(integer())
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  values <- runs$values
  left <- c(-Inf, utils::head(values, -1L))
  right <- c(utils::tail(values, -1L), -Inf)
  keep <- values > 0 & values > left & values > right
  as.integer(floor((starts[keep] + ends[keep]) / 2))
}
