#' Calculate mean squared deviation and its components
#'
#' @description
#' Compare paired predictions and observations using mean squared deviation
#' (MSD) and its decomposition into squared bias, nonunity slope, and lack of
#' correlation, following \insertCite{Gauch2003;textual}{rcaiman}.
#'
#' @details
#' Let \eqn{X} denote predictions and \eqn{Y} observations. MSD is
#' \eqn{\mathrm{mean}((X-Y)^2)}. Its decomposition is
#' \deqn{\mathrm{MSD} = \mathrm{SB} + \mathrm{NU} + \mathrm{LC},}
#' with
#' \deqn{\mathrm{SB} = (\bar{X}-\bar{Y})^2,}
#' \deqn{\mathrm{NU} = (1-b)^2\,\mathrm{mean}((X-\bar{X})^2),}
#' \deqn{\mathrm{LC} = \mathrm{mean}((Y-\widehat{Y})^2).}
#' Here \eqn{b} is the slope of observations regressed on predictions, with an
#' intercept, and \eqn{\widehat{Y}} is the fitted value of that regression.
#' The calculation uses centered moments with divisor \eqn{n}, rather than
#' sample variances with divisor \eqn{n-1}.
#'
#' SB measures the squared difference between means. NU measures departure
#' of the regression slope from one. LC measures the remaining variation of
#' observations around the regression line. When correlation is defined,
#' LC also equals \eqn{(1-r^2)\,\mathrm{mean}((Y-\bar{Y})^2)}.
#' All components have the squared units of the supplied values.
#'
#' Swapping predictions and observations preserves MSD and SB, but can change
#' NU and LC. The regression orientation used here agrees with
#' \insertCite{Pineiro2008;textual}{rcaiman}, who recommend regressing observations
#' on predictions when evaluating slope and intercept against the 1:1 line.
#'
#' If predictions are constant, their regression slope is undefined;
#' the decomposition is extended by setting NU to zero and LC to
#' \eqn{\mathrm{mean}((Y-\bar{Y})^2)}. With one pair, MSD equals SB and both
#' other components are zero. Equality between MSD and the sum of its
#' components holds up to floating-point rounding.
#'
#' The function summarizes the supplied pairs. Independent model validation
#' depends on how those predictions and observations were obtained.
#'
#' @param predicted,observed nonempty finite numeric vectors of equal length.
#'   Values are paired by position and must use the same units. No recycling
#'   or removal of missing values is performed.
#'
#' @return Numeric vector with named elements `MSD`, `SB` (squared bias),
#'   `NU` (nonunity slope), and `LC` (lack of correlation).
#'
#' @references \insertAllCited{}
#' @seealso [fit_cie_model()], [validate_cie_model()]
#' @export
#' @examples
#' predicted <- c(2, 4, 7, 8)
#' observed <- c(1, 3, 3, 5)
#' comparison <- calc_msd(predicted, observed)
#' comparison
#' comparison["MSD"]
#' sum(comparison[c("SB", "NU", "LC")])
#' calc_msd(observed, predicted)
calc_msd <- function(predicted, observed) {
  .check_vector(predicted, "numeric")
  .check_vector(observed, "numeric")
  if (!is.null(dim(predicted)) || !is.null(dim(observed)) ||
      !length(predicted) || length(predicted) != length(observed)) {
    stop("`predicted` and `observed` must be nonempty vectors of equal length.",
         call. = FALSE)
  }
  x <- predicted - mean(predicted)
  y <- observed - mean(observed)
  vx <- mean(x^2)
  slope <- if (vx == 0) 0 else mean(x * y) / vx
  c(MSD = mean((predicted - observed)^2),
    SB = (mean(predicted) - mean(observed))^2,
    NU = (1 - slope)^2 * vx,
    LC = mean((y - slope * x)^2))
}
