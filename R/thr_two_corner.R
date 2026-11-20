#' Compute thresholds with the two-corner method
#'
#' @description
#' Compute thresholds between two histogram corners using the method of
#' \insertCite{Macfarlane2011;textual}{rcaiman}, based on
#' \insertCite{Rosin2001;textual}{rcaiman}. Histogram preparation and peak
#' selection are separate steps.
#'
#' @details
#' The lower peak is used to find a corner toward higher values; the upper
#' peak is used to find a corner toward lower values. Each construction follows
#' [thr_rosin()]. With `slope_reduction = TRUE`, a peak above the mean bin
#' frequency is lowered to that mean when constructing its reference line,
#' following \insertCite{Macfarlane2011;textual}{rcaiman}.
#'
#' All thresholds are computed from `h`, including the minimum between the
#' corners (`ts`) \insertCite{Macfarlane2011}{rcaiman}. To use a smoothed curve,
#' pass the result of [smooth_histogram()] as `h`. Smoothing is a separate
#' choice made by the researcher; it is not applied internally. Peaks may be
#' located on another version of the histogram, provided the bin positions are
#' unchanged. If there is no bin strictly between the corners, `ts` is `NA` with
#' a warning; the other thresholds are still returned.
#'
#' [find_macfarlane_peaks()] requires exactly 256 bins because its search
#' windows were defined on that scale. Prepare them explicitly with
#' `build_histogram(x, n_bins = 256)`. This is independent of image bit depth.
#' When peaks come from another procedure, the corner calculation itself
#' accepts other numbers of equally spaced bins. Returned values always use
#' the original observation units.
#'
#' @section Exploring modality and locating peaks:
#' [find_histogram_peaks()] identifies local maxima in the supplied histogram.
#' Comparing their number and positions after [smooth_histogram()] with
#' different values of `sigma` helps distinguish persistent features from
#' small fluctuations. This describes the curve; it is not a statistical test
#' of the number of modes in the underlying population. Weighted histograms
#' can be explored in exactly the same way.
#'
#' For unweighted observations, [multimode::modetest()] offers formal tests.
#' With `mod0 = 1`, it tests unimodality against more than one mode.
#'
#' [multimode::locmodes()] estimates locations for a specified number of modes;
#' `mod0 = 2` requests two rather than testing for them. Its `locations` vector
#' alternates modes and antimodes: positions 1 and 3 are the two modes. Map
#' these locations to distinct histogram bin indices before supplying `peaks`,
#' as shown below. Neither function is called internally by `thr_two_corner()`.
#'
#' These `multimode` functions can receive the original unweighted observations,
#' not histogram frequencies. The use of solid-angle weights are not supported
#' by that alternative workflow.
#'
#' @inheritParams thr_rosin
#' @param peaks integer-like numeric vector of length two. Distinct bin indices
#'   in increasing order, starting at 1. They must have positive frequency in
#'   the curve used for corner detection. They need not be local maxima of `h`:
#'   for example, they may have been located before smoothing.
#' @param slope_reduction logical vector of length one. Apply the modification
#'   described in *Details*. Default `FALSE`.
#'
#' @return A list with `lp` and `up` (peak locations), `lc` and `uc` (corner
#'   locations), thresholds `tl`, `tm`, and `th` at 25%, 50%, and 75% of the
#'   interval between the corners, and `ts` (the first minimum-frequency bin
#'   strictly between the corners). All values use the original observation
#'   units. An error is raised if a corner cannot be found or `lc > uc`.
#'
#' @seealso [build_histogram()], [smooth_histogram()],
#'   [find_macfarlane_peaks()], [find_histogram_peaks()], [thr_rosin()],
#'   [binarize_by_region()], [apply_by_direction()]
#' @references \insertAllCited{}
#' @export
#' @examples
#' x <- c(qnorm(ppoints(2000), 30, 7), qnorm(ppoints(1800), 210, 9),
#'        seq(0, 255, length.out = 200))
#' h <- build_histogram(x, n_bins = 256)
#' peaks <- find_macfarlane_peaks(h)
#' h_smooth <- smooth_histogram(h, sigma = 2)
#' thr_two_corner(h_smooth, peaks, slope_reduction = TRUE, diagnose = TRUE)
#'
#' local_peaks <- find_histogram_peaks(h_smooth)
#' h$mids[local_peaks]
#'
#' # An alternative using the original, unweighted observations:
#' \dontrun{
#' assessment <- multimode::modetest(x, mod0 = 1)
#' assessment
#' modes <- multimode::locmodes(x, mod0 = 2, display = TRUE)
#' locations <- modes$locations[c(1, 3)]
#' peaks <- vapply(locations, function(value) {
#'   which.min(abs(h$mids - value))
#' }, integer(1))
#' thr_two_corner(h_smooth, peaks)
#' }
thr_two_corner <- function(h, peaks, slope_reduction = FALSE,
                           diagnose = FALSE, adjust_par = TRUE) {
  .check_threshold_histogram(h)
  .check_histogram_peaks(peaks, h, 2L)
  .check_vector(slope_reduction, "logical", 1)
  .check_vector(diagnose, "logical", 1)
  .check_vector(adjust_par, "logical", 1)
  left <- .histogram_corner(h, peaks[1], "right", slope_reduction)
  right <- .histogram_corner(h, peaks[2], "left", slope_reduction)
  lc <- unname(h$mids[left$index])
  uc <- unname(h$mids[right$index])
  if (lc > uc) stop("The lower corner is higher than the upper corner.", call. = FALSE)
  between <- which(seq_along(h$counts) > left$index & seq_along(h$counts) < right$index)
  ts <- NA_real_
  if (length(between)) {
    ts <- unname(h$mids[between[which.min(h$counts[between])]])
  } else {
    warning("No bin lies strictly between the corners; `ts` is NA.", call. = FALSE)
  }
  if (diagnose) {
    if (adjust_par) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar))
      graphics::par(mar = c(4, 4, 2, 1))
    }
    .plot_histogram_corners(h, list(left, right), c("right", "left"),
                           "Two-corner thresholds")
  }
  list(lp = unname(h$mids[peaks[1]]), lc = lc,
       tl = lc + (uc - lc) * .25, tm = lc + (uc - lc) * .5,
       th = lc + (uc - lc) * .75, uc = uc,
       up = unname(h$mids[peaks[2]]), ts = ts)
}
