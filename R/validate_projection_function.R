#' Validate radial polynomial projection function
#'
#' @description
#' Validation refers to compliance with the normalization constraints
#' required by the package for radial polynomial projection models.
#' It does not assess physical calibration correctness.
#'
#' @details
#' The public interface evaluates the model at 0 and 90 deg,
#' but the polynomial itself is defined in radians. Internally,
#' angles are converted to radians before evaluating the projection
#' function (see [calc_relative_radius()]).
#'
#' The validation distinguishes between structural validity and
#' operational safety:
#'
#' \enumerate{
#' \item Structural validity:
#' \itemize{
#' \item The projection model must satisfy the theoretical boundary conditions:
#' relative radius (R) at 0 deg equal to cero and R at 90 deg equal to one. Large
#' deviations (greater than `tol`) indicate an invalid calibration and require
#' recalibration.
#' \item The projection function must be monotonically increasing over the
#' interval `[0, 90 deg]`, i.e., the relative radius must not decrease as the
#' zenith angle increases. Violations of this condition indicate a geometrically
#' inconsistent model (e.g., radial folding), where multiple zenith angles map
#' to the same radius. Such models are considered structurally invalid and
#' require recalibration.
#' }
#' \item Operational safety: If R(90 deg) < 1 but the deviation is within
#' `tol`, the model is considered structurally valid but operationally unsafe
#' for segmentation algorithms that do not contemplate values above 90 deg. In
#' this case, [fix_projection_coef()] is recommended.
#' }
#'
#' @param projection_coef numeric vector. Polynomial coefficients defining
#'   \eqn{r(\theta) = \sum a_p \theta^p}, where \eqn{\theta} is expressed in
#'   radians.
#' @param tol numeric vector of length one. Maximum absolute deviation allowed
#'   at the boundaries before the model is considered structurally invalid.
#'
#' @return
#' A named list with elements:
#' \describe{
#'   \item{status}{Character string taking one of
#'     `"valid"`, `"structural_invalid"`, or `"underestimation"`.}
#'   \item{val_0}{Value of the projection function at 0 deg.}
#'   \item{val_90}{Value of the projection function at 90 deg.}
#'   \item{err_90}{Deviation from 1 at 90 deg.}
#' }
#'
#' @export
#'
#' @examples
#' validate_projection_function(c(0.6890,  0.0131, -0.0297), tol = 1e-3)
#' validate_projection_function(c(0.6890,  0.0131, -0.02), tol = 1e-3)
validate_projection_function <- function(projection_coef, tol = 1e-3) {

    .check_vector(tol, "numeric", 1, sign = "positive")
    if (tol > 1e-2) stop( "`tol` cannot be greather than 0.01.")

    val_0  <- calc_relative_radius(0,  projection_coef)
    val_90 <- calc_relative_radius(90, projection_coef)
    err_90 <- val_90 - 1

    # Structural invalid at 0 deg
    if (abs(val_0) > tol) {
      message(
        sprintf(
          paste0(
            "Projection function is structurally invalid at 0 deg.\n",
            "  f(0) = %.10f\n",
            "  Allowed deviation = +/- %.10f.\n",
            "Recalibration is required."
          ),
          val_0,
          tol
        )
      )
      return(list(
        status = "structural_invalid",
        val_0  = val_0,
        val_90 = val_90,
        err_90 = err_90
      ))
    }

    # Structural invalid at 90 deg
    if (abs(err_90) > tol) {
      message(
        sprintf(
          paste0(
            "Projection function is structurally invalid at 90 deg.\n",
            "  f(90 deg) = %.10f\n",
            "  Allowed deviation = +/- %.10f.\n",
            "Recalibration is required."
          ),
          val_90,
          tol
        )
      )
      return(list(
        status = "structural_invalid",
        val_0  = val_0,
        val_90 = val_90,
        err_90 = err_90
      ))
    }

    # Underestimation (operational issue only)
    if (err_90 < 0) {
      message(
        sprintf(
          paste0(
            "Projection function underestimates relative radius at 90 deg.\n",
            "  f(90 deg) = %.10f\n",
            "This may trigger hard-to-trace errors in downstream functions.\n",
            "Consider using 'fix_projection_coef()'."
          ),
          val_90
        )
      )
      return(list(
        status = "underestimation",
        val_0  = val_0,
        val_90 = val_90,
        err_90 = err_90
      ))
    }

    # Monotonicity check
    theta_seq <- seq(0, 90, length.out = 1000)
    r_vals <- calc_relative_radius(theta_seq, projection_coef)

    if (any(diff(r_vals) < 0)) {
      message(
        paste0(
          "Projection function is not monotonic increasing.\n",
          "This indicates a structurally invalid geometric model.\n",
          "Recalibration is required."
        )
      )
      return(list(
        status = "structural_invalid",
        val_0  = val_0,
        val_90 = val_90,
        err_90 = err_90
      ))
    }

    # Fully valid
    list(
      status = "valid",
      val_0  = val_0,
      val_90 = val_90,
      err_90 = err_90
    )
}
