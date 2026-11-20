#' Find histogram peaks using Macfarlane's windows
#'
#' @description
#' Locate two peaks using the windowed search of
#' \insertCite{Macfarlane2011;textual}{rcaiman}, separately from the
#' two-corner threshold construction.
#'
#' @details
#' This search requires exactly 256 equally spaced bins. It uses positions
#' 0--255 internally, regardless of the original data units. Initial windows
#' span positions 5--55 and 200--250. The lower window expands upward and the
#' upper window downward in steps of 25 until each maximum is at least ten
#' positions from the respective inner window limit. Ties select the lower
#' position. The peaks must be at least 20 positions apart and have positive
#' frequency. A failed search raises an error; it is not a test of unimodality.
#'
#' Use [build_histogram()] with `n_bins = 256` and supply the original,
#' unsmoothed histogram for the original Macfarlane workflow. Howerver, other
#' workflow are supported such as using sums of weights instead of direct counts
#' to compute frequencies.
#'
#' @inheritParams smooth_histogram
#'
#' @return Integer vector of length two, containing the lower and upper peak
#'   bin indices, starting at 1. Use `h$mids[peaks]` to obtain their values.
#'
#' @seealso [thr_two_corner()], [find_histogram_peaks()]
#'
#' @references \insertAllCited{}
#'
#' @export
#' @examples
#' x <- c(qnorm(ppoints(2000), 30, 7), qnorm(ppoints(1800), 210, 9),
#'        seq(0, 255, length.out = 200))
#' h <- build_histogram(x, n_bins = 256)
#' peaks <- find_macfarlane_peaks(h)
#' h$mids[peaks]
#' thr_two_corner(h, peaks, slope_reduction = TRUE)
find_macfarlane_peaks <- function(h) {
  .check_threshold_histogram(h)
  if (length(h$counts) != 256L) {
    stop("Macfarlane peak detection requires 256 bins; prepare the histogram with n_bins = 256.",
         call. = FALSE)
  }
  search <- function(lower, upper, side) {
    repeat {
      positions <- seq.int(lower, upper)
      peak <- positions[which.max(h$counts[positions + 1L])]
      margin <- if (side == "lower") upper - peak else peak - lower
      if (h$counts[peak + 1L] > 0 && margin >= 10) return(peak + 1L)
      if (side == "lower") {
        if (upper == 255) break
        upper <- min(255, upper + 25L)
      } else {
        if (lower == 0) break
        lower <- max(0, lower - 25L)
      }
    }
    stop("Macfarlane's ", side, " peak window could not locate a suitable peak.", call. = FALSE)
  }
  peaks <- c(search(5L, 55L, "lower"), search(200L, 250L, "upper"))
  if (diff(peaks) < 20) {
    stop("Macfarlane's peak windows did not yield peaks at least 20 bins apart.", call. = FALSE)
  }
  as.integer(peaks)
}
