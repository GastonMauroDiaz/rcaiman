#' Fix lens projection coefficients
#'
#' @description
#' Adjust the last coefficient of a radial polynomial projection model to
#' enforce the normalization constraint at 90 degrees (i.e., r(90°) >= 1). The
#' adjustment addresses operational boundary compliance required by the package.
#' It does not correct physical calibration errors.
#'
#' @details
#' This function applies exclusively to radial projection models of the form
#' \eqn{r(\theta) = \sum a_p \theta^p}, where \eqn{\theta} is expressed in
#' radians.
#'
#' This function relies on [validate_projection_function()] to classify
#' the projection model:
#'
#' \itemize{
#'  \item If the model is structurally invalid (large boundary deviation),
#'    an error is raised and recalibration is required.
#'  \item If the model is structurally valid but predicts a value below 1
#'    at 90 deg, the last coefficient is increased incrementally
#'    until validation succeeds.
#' }
#'
#' The function determines the correction direction empirically by
#' probing the local effect of the last coefficient on the predicted
#' value at 90 deg. It assumes a locally monotonic relationship
#' between the last coefficient and the boundary value.
#'
#' @inheritParams validate_projection_function
#'
#' @param tol numeric vector of length one. Structural tolerance passed to
#'   [validate_projection_function()].
#'
#' @return
#' A numeric vector of projection coefficients. The input vector is
#' returned unchanged if no correction is required.
#'
#' @seealso [validate_projection_function()]
#'
#' @export
fix_projection_coef <- function(projection_coef, tol = 1e-3) {

  res <- suppressMessages(validate_projection_function(projection_coef, tol))

  if (res$status == "structural_invalid") {
    stop(
      sprintf(
        paste0(
          "Projection model is structurally invalid under tolerance %.1e.\n",
          "Recalibration is required."
        ),
        tol
      ),
      call. = FALSE
    )
  }

  if (res$status == "valid") {
    message("Nothing to fix")
    return(projection_coef)
  }

  # Only remaining case: underestimation
  message(
    paste0(
      "Projection function underestimates relative radius at 90 deg.\n",
      "This may trigger hard-to-trace errors in downstream functions.\n",
      "Returning fixed coeficients."
    )
  )

  last <- length(projection_coef)
  step <- tol / 100

  val_base <- calc_relative_radius(90, projection_coef)

  # Probe direction
  coef_test <- projection_coef
  coef_test[last] <- coef_test[last] + step
  val_test <- calc_relative_radius(90, coef_test)

  direction <- sign(val_test - val_base)

  if (direction == 0) {
    stop("Last coefficient does not affect value at 90 deg.")
  }
  max_iter <- 50L
  iter <- 0L

  repeat {

    iter <- iter + 1L

    projection_coef[last] <-
      projection_coef[last] + direction * step

    res <- suppressMessages(validate_projection_function(projection_coef, tol))

    if (res$status == "valid") break

    if (iter >= max_iter) {
      stop(
        paste0(
          "Failed to fix projection coefficients after ", max_iter, " iterations.\n",
          "Consider relaxing the tolerance (tol)."
        ),
        call. = FALSE
      )
    }
  }

  projection_coef
}
