#' Adjust the last radial projection coefficient
#'
#' @description
#' Correct a small deficit of relative radius at the normalization angle by
#' incrementing the last polynomial coefficient. Other coefficients stay fixed.
#'
#' @inheritParams validate_projection_function
#'
#' @details
#' [validate_projection_function()] must classify the input as valid or as an
#' underestimation within `tol`. Structurally invalid models are rejected with
#' the corresponding diagnostic; this function does not repair radial folding.
#'
#' The last coefficient is increased iteratively until the reference radius
#' reaches one, within floating-point precision. Each increment targets a change
#' of `tol / 100` in relative radius at `normalization_angle`, so the step remains
#' useful when that angle or the polynomial degree changes. The complete model
#' is validated after each increment, including invertibility up to
#' `max_off_axis_angle`. Failure to obtain a valid result raises an error.
#'
#' This preserves the existing correction mechanism for rounded coefficients.
#' It changes the curve through its highest-order term; it does not improve the
#' experimental fit. [normalize_projection_coef()] is a separate operation that
#' rescales all coefficients to change the normalization reference.
#'
#' @return Numeric vector of coefficients. A valid input is returned unchanged,
#'   including its names. An underestimated model is returned only after the
#'   corrected coefficients pass validation.
#' @seealso [validate_projection_function()]
#' @export
#' @examples
#' projection_coef <- c(0.702, 0.0434, -0.0543)
#' corrected <- fix_projection_coef(projection_coef)
#' rbind(original = projection_coef, corrected = corrected)
#' calc_relative_radius(90, corrected)
fix_projection_coef <- function(projection_coef, tol = 1e-3,
                                 normalization_angle = 90,
                                 max_off_axis_angle = 90) {
  validate <- function(x) validate_projection_function(
    x, tol, normalization_angle, max_off_axis_angle)
  result <- validate(projection_coef)
  if (result$status == "structural_invalid") stop(result$reason, call. = FALSE)
  if (result$status == "valid") return(projection_coef)

  last <- length(projection_coef)
  step <- (tol / 100) / (normalization_angle * pi / 180)^last
  if (!is.finite(step) || step <= 0) {
    stop("Cannot represent a finite correction step for the last coefficient.", call. = FALSE)
  }
  for (iteration in seq_len(200)) {
    previous <- projection_coef[last]
    projection_coef[last] <- previous + step
    if (!is.finite(projection_coef[last]) || projection_coef[last] == previous) {
      stop("The correction step cannot change the last coefficient at this precision.", call. = FALSE)
    }
    result <- validate(projection_coef)
    if (result$status == "valid") return(projection_coef)
    if (result$status == "structural_invalid") {
      stop("The last-coefficient adjustment did not preserve model validity: ",
           result$reason, call. = FALSE)
    }
  }
  stop("No valid correction was found after 200 increments.", call. = FALSE)
}
