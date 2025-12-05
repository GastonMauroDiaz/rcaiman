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
#' An observation \eqn{x_i} is considered outlying if it satisfies:
#'
#' \deqn{M - \mathrm{laxity} \times \mathrm{MAD} < x_i <
#'       M + \mathrm{laxity} \times \mathrm{MAD}}
#'
#' where \eqn{M} is the median of `x` and \eqn{MAD} is its median absolute
#' deviation. The `laxity` parameter multiplies MAD to set the classification
#' threshold.
#'
#' The `cutoff_side` argument controls which tail(s) are evaluated:
#' `"both"` (default) applies a two-sided test, `"left"` applies the left-tail
#' test only, and `"right"` applies the right-tail test only.
#'
#' If the MAD is zero (no dispersion), all values are flagged as non-outliers. `NA`
#' values in `x` are treated as non-outliers.
#'
#' @param x numeric vector. Values to evaluate.
#' @param y optinal numeric vector. Values from which median and MAD are
#'   calculated. If `NUll` (default) median and MAD will be calculated from `x`.
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
  mad_val <- stats::mad(y, constant = 1, na.rm = TRUE)

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
