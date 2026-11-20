#' Validate a radial projection model
#'
#' @description
#' Check normalization and invertibility of a radial projection model. These
#' checks assess mathematical consistency, not experimental calibration accuracy.
#'
#' @param projection_coef numeric vector of coefficients for the radial model
#'   described in [set_internal_geometry()]. The polynomial uses off-axis angles
#'   in radians and has no constant term.
#' @param tol positive numeric value, at most 0.01. Maximum absolute deviation
#'   from one at `normalization_angle`. Defaults to 0.001.
#' @param normalization_angle positive numeric value, at most 180. Off-axis
#'   angle in degrees at which the relative radius should equal one. Defaults
#'   to 90, independently of the inversion limit.
#' @param max_off_axis_angle positive numeric value, at most 180. Upper limit
#'   in degrees of the interval checked for invertibility. Defaults to 90.
#'
#' @details
#' The model must be finite and strictly increasing from zero to
#' `max_off_axis_angle`. The derivative is checked between its real roots,
#' using the same check as [calc_off_axis_angle()]. This check also applies
#' when the normalization radius is underestimated. The absent constant term
#' ensures a relative radius of zero on the optical axis.
#'
#' A reference radius within `tol` but below one is classified as
#' `"underestimation"`. When normalization and inversion angles coincide,
#' radii between that value and one have no inverse within the allowed interval.
#' Accepting an absolute normalization error alone would not prevent unexpected
#' missing directions near the nominal circle boundary. [fix_projection_coef()]
#' adjusts the last coefficient to resolve this small deficit.
#'
#' Differences within eight machine epsilons relative to the reference radius
#' are treated as floating-point roundoff, consistently with the radial inverse.
#' This is distinct from `tol`, which accommodates rounded model coefficients.
#' [set_internal_geometry()] uses this validation with `tol = 0.01` and rejects
#' both invalid and underestimated models without changing their coefficients.
#' A calculation limit below the normalization angle still limits coverage:
#' validation does not extend the requested angular domain.
#'
#' @return Named list with `status` (`"valid"`, `"structural_invalid"`, or
#'   `"underestimation"`), `val_0`, `val_reference` (radius at
#'   `normalization_angle`), `err_reference` (that radius minus one), and
#'   `reason` (a diagnostic string, or `NULL` for a valid model). Malformed
#'   arguments raise errors; model failures are reported through this list.
#' @seealso [fix_projection_coef()], [normalize_projection_coef()]
#' @export
#' @examples
#' # Rounded coefficients from a lens calibration.
#' projection_coef <- c(0.702, 0.0434, -0.0543)
#' validate_projection_function(projection_coef)
#' corrected <- fix_projection_coef(projection_coef)
#' validate_projection_function(corrected)
#'
#' # Normalization and calculation limit need not coincide.
#' validate_projection_function(3 / pi, normalization_angle = 60,
#'                              max_off_axis_angle = 120)
validate_projection_function <- function(projection_coef, tol = 1e-3,
                                          normalization_angle = 90,
                                          max_off_axis_angle = 90) {
  .check_vector(projection_coef, "numeric")
  if (!length(projection_coef)) stop("`projection_coef` cannot be empty.", call. = FALSE)
  .check_vector(tol, "numeric", 1, sign = "positive")
  if (tol > 0.01) stop("`tol` must not exceed 0.01.", call. = FALSE)
  for (name in c("normalization_angle", "max_off_axis_angle")) {
    value <- get(name)
    .check_vector(value, "numeric", 1, sign = "positive", name = name)
    if (value > 180) stop("`", name, "` must not exceed 180 degrees.", call. = FALSE)
  }
  reference <- unname(.calc_relative_radius(normalization_angle * pi / 180, projection_coef))
  error <- reference - 1
  inversion_error <- tryCatch({
    .check_radial_inversion(projection_coef, max_off_axis_angle)
    NULL
  }, error = function(e) conditionMessage(e))
  status <- "valid"
  reason <- NULL
  if (!is.null(inversion_error)) {
    status <- "structural_invalid"
    reason <- inversion_error
  } else if (!is.finite(reference) || abs(error) > tol) {
    status <- "structural_invalid"
    reason <- paste0("The radial polynomial must equal one at `normalization_angle` ",
                     "within tolerance ", tol, "; obtained ", format(reference, digits = 12), ".")
  } else if (error < -8 * .Machine$double.eps * abs(reference)) {
    status <- "underestimation"
    reason <- paste0("The relative radius at `normalization_angle` is below one (",
                     format(reference, digits = 12), "). Use `fix_projection_coef()` ",
                     "with tol = ", tol, " and the same normalization and calculation angles.")
  }
  list(status = status, val_0 = 0, val_reference = reference,
       err_reference = error, reason = reason)
}
