.calc_relative_radius <- function(angle, projection_coef) {
  # Horner method
  n <- length(projection_coef)

  R <- rep(0, length(angle))
  for (i in n:1) {
    R <- (R + projection_coef[i]) * angle
  }

  R
}

#' Calculate relative radius (radial projection model)
#'
#' @description
#' Convert zenith angles to normalized radial distance under a
#' radial projection model with circular symmetry.
#'
#' @details
#' This function applies to hemispherical imaging systems where the mapping
#' between zenith angle and image position is radially symmetric, such as in
#' fisheye lenses.
#'
#' The projection function is defined as:
#' \deqn{
#' r_{rel}(\theta) = \sum_{p=1}^{n} a_p \theta^p
#' }
#' where \eqn{r_{rel}} is the relative radius (radius expresed as a fraction of
#' the radius at 90 deg), \eqn{\theta} is the zenith angle in radians and
#' `projection_coef = c(a_1, a_2, ..., a_n)`. The element at position `p`
#' multiplies \eqn{\theta^p}.
#'
#' This function can manage `angle` units to reflects common practice in
#' hemispherical photography software, where fitting is performed in radians but
#' operational workflows use degrees.
#'
#' Under a valid calibration (see [validate_projection_function()]), values are
#' expected to lie in \[0, 1\], with 0 at the zenith (0 deg) and 1 at the
#' horizon (90 deg).
#'
#' @param angle numeric vector. Zenith angles expressed as declared in `units`.
#' @param projection_coef numeric vector. Polynomial coefficients defining
#'   the radial projection model evaluated in radians.
#' @param units character vector of length one. Either `"radian"` or `"degree"`,
#'   for declaring the units in which `angle` is provided.
#'
#' @return Numeric vector of the same length as `angle`.
#'
#' @note
#' If it is forced to extrapolate by inputing values above 90 deg, the
#' [equidistant_coef] is applied to guarantee monotonicity.
#'
#' @export
#'
#' @examples
#' calc_relative_radius(45, equidistant_coef)
#' calc_relative_radius(45, equal_area_coef)
#' calc_relative_radius(45, orthographic_coef)
calc_relative_radius <- function(angle, projection_coef, units = "degree") {
  .assert_choice(units, c("radian", "degree"))
  .check_vector(angle, "numeric", sign = "any")
  .check_vector(projection_coef, "numeric", sign = "any")

  if (units == "degree") {
    angle <- .degree2radian(angle)
  }

  i <- angle > pi/2
  relative_radiance <- if (any(i)) {
    powers <- equidistant_coef #to extrapolate
    c(.calc_relative_radius(angle[!i], projection_coef),
      .calc_relative_radius(angle[i], equidistant_coef))
  } else {
    .calc_relative_radius(angle, projection_coef)
  }
  relative_radiance
}
