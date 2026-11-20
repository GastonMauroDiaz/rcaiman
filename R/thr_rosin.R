#' Compute Rosin's unimodal threshold
#'
#' @description
#' Find a geometric corner on a selected tail of a histogram using
#' \insertCite{Rosin2001;textual}{rcaiman}. The peak is supplied separately,
#' allowing the researcher to choose how to locate it and prepare the histogram.
#'
#' @details
#' A line joins the supplied peak to a zero-frequency position one bin beyond
#' the histogram on the selected side. The threshold is the bin below that
#' line with the greatest perpendicular distance to it. Tied corners use the
#' first bin encountered from the peak toward that end. For a left tail, the
#' histogram is reflected before finding the corner.
#'
#' No smoothing or peak search is performed internally. Use [build_histogram()]
#' to prepare observations, optionally with weights, and [smooth_histogram()]
#' to smooth the resulting curve. A dominant peak can be selected with
#' `which.max(h$density)`; [find_histogram_peaks()] allows inspection of local
#' maxima.
#'
#' Unlike [find_macfarlane_peaks()], this calculation accepts any number of
#' equally spaced bins. Histogram resolution and smoothing affect the result.
#' Multiplying all frequencies by the same positive constant does not change it.
#'
#' Complementary gradients obtained with [complementary_gradients()] typically
#' have unimodal distributions. A dominant canopy population may leave a right
#' tail toward sky values. A dominant sky population may instead leave a left
#' tail toward canopy values. Choosing a tail does not change the gradient's
#' polarity.
#'
#' @param h histogram from [build_histogram()] or [smooth_histogram()], or a
#'   list with numeric `breaks`, `mids`, and nonnegative `counts`. Bins must be
#'   equally spaced, with matching midpoints and positive total frequency.
#' @param peak integer-like numeric vector of length one. Bin index of the
#'   peak, starting at 1, with positive frequency. This is an index, not the
#'   value at that position; `h$mids[peak]` gives the latter.
#' @param tail character vector of length one. `"right"` (default) searches
#'   above the peak; `"left"` searches below it.
#' @param diagnose logical vector of length one. If `TRUE`, plot the histogram,
#'   reference line, and corner. Horizontal tick labels and the returned
#'   threshold use the units of the observations.
#' @param adjust_par logical vector of length one. If `TRUE` (default), adjust
#'   plot margins for the diagnostic diagram and restore them afterward.
#'
#' @return Numeric vector of length one: the midpoint of the selected bin,
#'   in the original units of the observations. An error is raised if there
#'   is no corner on the selected side.
#'
#' @seealso [thr_two_corner()], [find_histogram_peaks()], [thr_isodata()],
#'   [apply_by_direction()], [complementary_gradients()], [binarize_with_thr()]
#' @references \insertAllCited{}
#' @export
#' @examples
#' x <- 2 * qbeta(ppoints(2000), shape1 = 2, shape2 = 10) - 1
#' h <- build_histogram(x)
#' h_smooth <- smooth_histogram(h)
#' peak <- which.max(h_smooth$counts)
#' thr_rosin(h_smooth, peak, diagnose = TRUE)
#'
thr_rosin <- function(h,
                      peak,
                      tail = "right",
                      diagnose = FALSE,
                      adjust_par = TRUE) {

  .check_threshold_histogram(h)
  .check_histogram_peaks(peak, h, 1L)
  .assert_choice(tail, c("right", "left"))
  .check_vector(diagnose, "logical", 1)
  corner <- .histogram_corner(h, peak, tail)
  threshold <- unname(h$mids[corner$index])
  if (diagnose) {
    if (adjust_par) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar))
      graphics::par(mar = c(4, 4, 2, 1))
    }
    .plot_histogram_corners(h, list(corner), tail,
      paste("Rosin threshold:", signif(threshold, 5)))
  }
  threshold
}
