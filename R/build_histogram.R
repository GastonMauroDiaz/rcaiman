#' Build a histogram with optional observation weights
#'
#' @description
#' Accumulate observations or their weights into an exact number of equally
#' spaced bins, retaining the units of the observations.
#'
#' @details
#' The breaks span the minimum and maximum retained values. Bins include their
#' lower boundary and exclude their upper boundary, except for the last bin,
#' which also includes the maximum. Values are not rounded to digital numbers.
#' `n_bins` is an exact count, not a suggestion to [graphics::hist()].
#'
#' Without weights, `counts` contains the number of observations per bin.
#' With weights, it contains their sums. Weights from [area_weight_image()]
#' therefore describe the distribution of solid angle rather than pixel counts.
#' Values of `x` are never multiplied or divided by weights. Zero weights and
#' pairs with missing values are removed before determining the range.
#'
#' The number of bins controls the resolution of the distribution, independently
#' of the bit depth of an image. The default 256 also meets the requirement of
#' [find_macfarlane_peaks()]. Changing the number of bins changes the width of
#' the Gaussian kernel used by [smooth_histogram()] for a fixed `sigma`.
#'
#' @inheritParams thr_isodata
#' @param n_bins integer-like numeric vector of length one, at least two.
#'   Exact number of bins. Default 256.
#'
#' @return A base R `histogram` object, with `breaks`, `mids`, `counts`,
#'   `density`, `xname`, and `equidist`.
#'   `density` integrates to one, including when counts are weighted. At least
#'   two distinct retained values and a finite representable range are required.
#'
#' @seealso [smooth_histogram()], [find_histogram_peaks()], [thr_rosin()],
#'   [thr_two_corner()]
#' @export
#' @examples
#' x <- c(0, 1, 2, 3)
#' h <- build_histogram(x, w = c(1, 1, 2, 4), n_bins = 3)
#' h$counts
#' plot(h, freq = FALSE)
build_histogram <- function(x, w = NULL, n_bins = 256) {
  xname <- paste(deparse(substitute(x)), collapse = " ")
  as_vector <- function(value, name) {
    if (is.matrix(value) || is.data.frame(value)) {
      if (ncol(value) != 1L) stop("`", name, "` must have a single column.", call. = FALSE)
      value <- if (is.data.frame(value)) value[[1]] else value[, 1]
    }
    .check_vector(value, "numeric", allow_na = TRUE, name = name)
    value
  }
  x <- as_vector(x, "x")
  if (is.null(w)) w <- rep(1, length(x)) else w <- as_vector(w, "w")
  .check_vector(w, "numeric", length(x), allow_na = TRUE,
                sign = "nonnegative", name = "w")
  .check_vector(n_bins, "numeric", 1)
  if (n_bins < 2 || n_bins != floor(n_bins) || n_bins >= .Machine$integer.max) {
    stop("`n_bins` must be an integer of at least two within R's integer range.", call. = FALSE)
  }
  keep <- !is.na(x) & !is.na(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) < 2L || min(x) == max(x)) {
    stop("At least two distinct retained values are required to build a histogram.", call. = FALSE)
  }
  span <- max(x) - min(x)
  if (!is.finite(span) || !is.finite(sum(w))) {
    stop("The data range and total weight must be finite.", call. = FALSE)
  }
  breaks <- seq(min(x), max(x), length.out = n_bins + 1L)
  if (any(diff(breaks) <= 0)) {
    stop("The data range cannot represent this many distinct bin boundaries.", call. = FALSE)
  }
  bins <- findInterval(x, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  sums <- rowsum(w, bins, reorder = FALSE)
  counts <- numeric(n_bins)
  counts[as.integer(rownames(sums))] <- sums[, 1]
  structure(list(breaks = breaks, counts = counts,
    density = counts / sum(counts) / diff(breaks),
    mids = utils::head(breaks, -1L) + diff(breaks) / 2,
    xname = xname, equidist = TRUE), class = "histogram")
}

.check_threshold_histogram <- function(h, name = deparse(substitute(h))) {
  if (!is.list(h) || !all(c("breaks", "mids", "counts") %in% names(h))) {
    stop("`", name, "` must contain histogram breaks, mids, and counts.", call. = FALSE)
  }
  n <- length(h$counts)
  .check_vector(h$counts, "numeric", sign = "nonnegative", name = paste0(name, "$counts"))
  .check_vector(h$breaks, "numeric", n + 1L, name = paste0(name, "$breaks"))
  .check_vector(h$mids, "numeric", n, name = paste0(name, "$mids"))
  widths <- diff(h$breaks)
  if (n < 2L || any(!is.finite(widths)) || any(widths <= 0) ||
      !is.finite(sum(h$counts)) || sum(h$counts) <= 0 ||
      any(abs(widths / mean(widths) - 1) > 1e-7) ||
      any(abs((h$mids - utils::head(h$breaks, -1L)) / widths - .5) > 1e-7)) {
    stop("`", name, "` must have equally spaced bins, matching midpoints, and positive total frequency.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.check_histogram_peaks <- function(peaks, h, n = length(peaks)) {
  .check_vector(peaks, "numeric", n, name = "peaks")
  if (!length(peaks) || any(peaks != floor(peaks)) ||
      any(peaks < 1 | peaks > length(h$counts)) || any(diff(peaks) <= 0)) {
    stop("Peaks must be distinct bin indices in increasing order, starting at 1.", call. = FALSE)
  }
  if (any(h$counts[peaks] <= 0)) stop("Peak bins must have positive frequency.", call. = FALSE)
  invisible(TRUE)
}
