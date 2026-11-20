#' Identify statistical outliers
#'
#' @description
#' Identify outlying values in a numeric vector using a robust MAD–median
#' rule and return a logical mask indicating which elements are outliers.
#'
#' @details
#' Apply the robust criterion from \insertCite{Leys2013;textual}{rcaiman} to
#' detect statistical outliers. The function returns a logical vector of the
#' same length as `x` with `TRUE` for detected outliers.
#'
#' For the two-sided test, an observation \eqn{x_i} is outlying when:
#'
#' \deqn{|x_i - M| > \mathrm{laxity} \times \mathrm{MAD},}
#'
#' where \eqn{M} is the median of the reference values `y` (or `x` when `y`
#' is `NULL`), and MAD is their median absolute deviation multiplied by
#' `1.4826`, as in [stats::mad()] by default. Values exactly on the threshold
#' are retained.
#'
#' The `cutoff_side` argument controls which tail(s) are evaluated:
#' `"both"` (default) applies a two-sided test, `"left"` applies the left-tail
#' test only, and `"right"` applies the right-tail test only.
#'
#' If MAD is zero or unavailable (an empty reference), all values are flagged
#' as non-outliers. Inputs must be finite; missing values are not accepted.
#'
#' @param x finite numeric vector. Values to evaluate.
#' @param y optional finite numeric vector. Reference values from which median
#'   and MAD are calculated. With `NULL` (default), use `x`.
#' @param laxity numeric vector of length one. Positive multiplier applied to
#'   MAD to set classification thresholds. Default is `2`.
#' @param cutoff_side character vector of length one. One of `"both"`,
#'   `"left"`, or `"right"`. Default is `"both"`.
#'
#' @return Logical vector of the same length as `x`. `TRUE` indicates the outliers.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' x <- c(rnorm(100), 10, 12)
#' i <- is_outlier(x, x, laxity = 2.5, cutoff_side = "both")
#' x[i]  # values flagged as outliers
#' x[!i] # filter outliers
is_outlier <- function(x, y = NULL, laxity = 2, cutoff_side = "both") {
  .check_vector(x, "numeric")
  .check_vector(y, "numeric", allow_null = TRUE)
  .check_vector(laxity, "numeric", 1, sign = "positive")
  .assert_choice(cutoff_side, c("both", "left", "right"))

  if (is.null(y)) y <- x

  med <- stats::median(y, na.rm = TRUE)
  mad_val <- stats::mad(y, na.rm = TRUE)

  # If dispersion is zero or NA, retain non-missing values (no evidence to remove)
  if (is.na(mad_val) || mad_val == 0) {
    return(rep(FALSE, length(x)))
  }

  deviation <- (x - med) / mad_val
  i <- switch(
    cutoff_side,
    right = deviation > laxity,
    left  = deviation < -laxity,
    both  = abs(deviation) > laxity
  )

  # For any NA in the comparison (e.g., due to NA in x), treat as retained
  i[is.na(i)] <- TRUE

  i
}
