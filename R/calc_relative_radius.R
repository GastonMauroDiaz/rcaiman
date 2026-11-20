.calc_relative_radius <- function(angle, projection_coef) {
  # Horner method
  R <- rep(0, length(angle))
  for (i in rev(seq_along(projection_coef))) {
    R <- (R + projection_coef[i]) * angle
  }
  R
}

#' Evaluate a radial projection in either direction
#'
#' @description
#' Convert between off-axis angles and relative radii using a radial
#' projection polynomial.
#'
#' @details
#' [calc_relative_radius()] evaluates the polynomial from angle to relative
#' radius. [calc_off_axis_angle()] numerically inverts the same polynomial
#' from relative radius to angle. The radial model and its angular convention
#' are described in [set_internal_geometry()].
#'
#' Relative radius is measured as a fraction of `normalization_radius`, the
#' radius in pixels at the chosen `normalization_angle`. That reference is
#' expressed by the coefficients: the polynomial equals one at the reference
#' angle. Use [normalize_projection_coef()] to change the reference.
#' Relative radii greater than one correspond to angles beyond that reference;
#' they do not necessarily lie beyond the field of view.
#'
#' To convert a measured pixel radius to an angle, first divide it by
#' `normalization_radius`. To convert an angle to a pixel radius, multiply
#' the relative radius returned by `calc_relative_radius()` by that same
#' reference radius.
#'
#' ### Evaluation and inversion
#'
#' The forward evaluator uses the supplied polynomial throughout, including
#' above 90 deg. It can also be used to inspect an unconstrained fit; evaluating
#' a polynomial does not establish that it is a valid lens projection.
#'
#' Inversion requires a strictly increasing polynomial on
#' `[0, max_off_axis_angle]`. This is the same condition and numerical inversion
#' used by [set_internal_geometry()]. Radii beyond the value of the polynomial
#' at the upper limit return `NA`. The limit is independent of normalization
#' and does not describe the range over which the lens was calibrated.
#'
#' @param angle numeric vector of off-axis angles, in `[0, 180]` degrees or
#'   `[0, pi]` radians, as declared by `units`. Missing values are retained.
#' @param projection_coef nonempty finite numeric vector. Coefficients in
#'   increasing order of degree, without an intercept. The polynomial is
#'   evaluated in radians, irrespective of `units`.
#' @param units character vector of length one, either `"degree"` or
#'   `"radian"`. Declares the input angle units for `calc_relative_radius()`
#'   and the output angle units for `calc_off_axis_angle()`.
#' @param relative_radius numeric vector of non-negative relative radii.
#'   Missing values are retained.
#' @param max_off_axis_angle positive numeric vector of length one, at most
#'   180. Upper limit for inversion, always in degrees. Defaults to 90.
#'
#' @return numeric vector with the same length, order, and names as the input.
#'   `calc_relative_radius()` returns dimensionless relative radii.
#'   `calc_off_axis_angle()` returns off-axis angles in the requested units,
#'   with `NA` for radii outside the inversion interval.
#'
#' @seealso [normalize_projection_coef()], [set_internal_geometry()]
#' @export
#' @examples
#' # A standard projection, normalized at 90 degrees.
#' rho <- calc_relative_radius(c(0, 30, 60, 90), equal_area_coef)
#' calc_off_axis_angle(rho, equal_area_coef)
#'
#' # Normalize at 60 degrees and calculate angles up to 120 degrees.
#' # This quadratic example illustrates the calculation, not a lens calibration.
#' projection_coef <- normalize_projection_coef(c(0, 1), 60)
#' rho <- calc_relative_radius(c(30, 60, 90, 120), projection_coef)
#' calc_off_axis_angle(rho, projection_coef, max_off_axis_angle = 120)
#'
#' # Estimate the direction corresponding to a measured pixel radius.
#' normalization_radius <- 300
#' measured_radius <- 400
#' calc_off_axis_angle(measured_radius / normalization_radius,
#'                     projection_coef, max_off_axis_angle = 120)
calc_relative_radius <- function(angle, projection_coef, units = "degree") {
  .assert_choice(units, c("radian", "degree"))
  .check_vector(angle, "numeric", sign = "nonnegative", allow_na = TRUE)
  .check_vector(projection_coef, "numeric")
  if (!length(projection_coef)) {
    stop("`projection_coef` cannot be empty.", call. = FALSE)
  }
  limit <- if (units == "degree") 180 else pi
  if (any(angle > limit, na.rm = TRUE)) {
    stop("Off-axis angles must not exceed 180 degrees (pi radians).",
         call. = FALSE)
  }
  if (units == "degree") angle <- .degree2radian(angle)
  result <- .calc_relative_radius(angle, projection_coef)
  names(result) <- names(angle)
  result
}

#' @rdname calc_relative_radius
#' @export
calc_off_axis_angle <- function(relative_radius, projection_coef,
                                max_off_axis_angle = 90, units = "degree") {
  .assert_choice(units, c("radian", "degree"))
  .check_vector(relative_radius, "numeric", sign = "nonnegative", allow_na = TRUE)
  .check_radial_inversion(projection_coef, max_off_axis_angle)
  angle <- .calc_off_axis_angle(relative_radius, projection_coef,
                                 max_off_axis_angle)
  if (units == "degree") angle <- angle * 180 / pi
  angle
}

.calc_off_axis_angle <- function(relative_radius, projection_coef, max_angle) {
  upper <- max_angle * pi / 180
  boundary <- .calc_relative_radius(upper, projection_coef)
  # Allow floating-point roundoff at the projected angular limit.
  at_boundary <- !is.na(relative_radius) &
    abs(relative_radius - boundary) <= 8 * .Machine$double.eps * boundary
  inside <- !is.na(relative_radius) &
    (relative_radius <= boundary | at_boundary)
  target <- relative_radius[inside]
  angle <- rep(NA_real_, length(relative_radius))
  if (length(target)) {
    low <- numeric(length(target))
    high <- rep(upper, length(target))
    # Vectorized bisection keeps inversion independent of image resolution.
    for (i in seq_len(40)) {
      middle <- (low + high) / 2
      below <- .calc_relative_radius(middle, projection_coef) < target
      low[below] <- middle[below]
      high[!below] <- middle[!below]
    }
    angle[inside] <- (low + high) / 2
    angle[which(relative_radius == 0)] <- 0
    angle[at_boundary] <- upper
  }
  names(angle) <- names(relative_radius)
  angle
}
