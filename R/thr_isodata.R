#' Compute IsoData threshold
#'
#' @description
#' Compute a threshold using the IsoData algorithm
#' \insertCite{isodata;textual}{rcaiman}, recommended by
#' \insertCite{Jonckheere2005;textual}{rcaiman}, with optional observation weights.
#'
#' @details
#' The unweighted calculation follows the IsoData method by Gabriel Landini,
#' as implemented in [autothresholdr::auto_thresh()], while accepting numeric
#' data over an arbitrary range. The weighted extension uses the same iterative
#' partition with weighted means.
#'
#' The initial threshold is the weighted mean of all retained values. At each
#' iteration, values at or below the threshold form one group and values above
#' it form the other. The new threshold is the midpoint of the two group means:
#' \deqn{T_{new} = (\mu_{low} + \mu_{high}) / 2,}
#' where each mean is computed using the weights of its own observations.
#' The midpoint is not weighted by the total weight of either group. Iteration
#' ends when the threshold stops changing. All retained observations are used;
#' no random subsampling is performed.
#'
#' @section Solid-angle weighting:
#' Weights change the contribution of each observation to the DN distribution,
#' not its DN. With weights from [area_weight_image()], each pixel contributes
#' according to the solid angle it represents. Dividing or multiplying image
#' values by those weights would instead change the radiometric data.
#'
#' Relative weights suffice: multiplying all weights by the same positive
#' constant leaves the threshold unchanged. Equal positive weights give the
#' unweighted result. The threshold remains in the units of `x` and can be
#' applied directly to the original values with [binarize_with_thr()]. For
#' angular windows, [apply_by_direction()] passes matching values and weights
#' directly to this function.
#'
#' @param x numeric vector, or a single-column numeric `matrix` or `data.frame`.
#'   Missing values are ignored; infinite values are not allowed.
#' @param w optional numeric vector, or a single-column numeric `matrix` or
#'   `data.frame`, with one weight per value in `x`. Weights must be nonnegative
#'   and finite, or `NA`. A missing value in either input excludes that pair;
#'   zero weights also exclude observations. With `NULL`, every observation
#'   has weight one. At least two observations with positive weights and
#'   nonmissing values must remain.
#'
#' @return Numeric vector of length one, in the units of `x`. If all retained
#'   values are equal, that value is returned with a warning.
#' @seealso [area_weight_image()], [apply_by_direction()], [binarize_with_thr()]
#' @export
#' @references \insertAllCited{}
#' @examples
#' # Weighted observations are not modified digital numbers.
#' x <- c(0, 2, 10, 12)
#' thr_isodata(x)
#' thr_isodata(x, w = c(1, 3, 1, 1))
#'
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, fill_na = "none")$Blue
#' r <- crop_caim(r, c(600, 800), 200, 200)
#'
#' thr <- thr_isodata(r[])
#' plot(binarize_with_thr(r, thr))
#' # See apply_by_direction() for thresholds in angular windows.
thr_isodata <- function(x, w = NULL) {
  as_vector <- function(value, name) {
    if (is.matrix(value) || is.data.frame(value)) {
      if (ncol(value) != 1L) {
        stop("`", name, "` must have a single column.", call. = FALSE)
      }
      value <- if (is.data.frame(value)) value[[1]] else value[, 1]
    }
    .check_vector(value, "numeric", allow_na = TRUE, name = name)
    value
  }
  x <- as_vector(x, "x")
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as_vector(w, "w")
    .check_vector(w, "numeric", length(x), allow_na = TRUE,
                  sign = "nonnegative", name = "w")
  }
  selected <- !is.na(x) & !is.na(w) & w > 0
  x <- x[selected]
  w <- w[selected]
  if (length(x) < 2L) {
    stop("At least two observations with positive weights and nonmissing values are required.",
         call. = FALSE)
  }
  if (all(x == x[1])) {
    warning("All retained values are equal; returning that value.", call. = FALSE)
    return(unname(x[1]))
  }
  # Only relative weights matter; rescale to avoid overflow in their sum.
  w <- w / max(w)
  threshold <- stats::weighted.mean(x, w)
  repeat {
    low <- x <= threshold
    mean_low <- stats::weighted.mean(x[low], w[low])
    mean_high <- stats::weighted.mean(x[!low], w[!low])
    updated <- mean_low / 2 + mean_high / 2
    if (!is.finite(updated)) {
      stop("Cannot compute a finite threshold from the supplied values and weights.",
           call. = FALSE)
    }
    if (updated == threshold) break
    threshold <- updated
  }
  unname(threshold)
}
