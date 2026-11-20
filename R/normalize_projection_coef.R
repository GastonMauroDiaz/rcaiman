#' Normalize radial projection coefficients
#'
#' @description
#' Rescale a radial projection polynomial to give a relative radius of one at
#' a chosen off-axis angle.
#'
#' @details
#' In the radial model of [set_internal_geometry()], a polynomial
#' \eqn{P(\theta') = \sum_{k=1}^{p} a_k (\theta')^k} predicts relative radius
#' from off-axis angle in radians. For a new reference angle
#' \eqn{\theta'_{norm}}, divide every coefficient by
#' \eqn{s = P(\theta'_{norm})}. The resulting polynomial is
#' \eqn{P_{new}(\theta') = P(\theta') / s} and equals one at that reference.
#' This changes the scale of relative radius while preserving the shape of
#' the projection curve.
#'
#' ### Standard projections
#'
#' [equidistant_coef], [equal_area_coef], and [orthographic_coef] are supplied
#' with a normalization reference of 90 deg. Use this function to express any
#' of them relative to a different angle. The equidistant model is exact;
#' the equal-area and orthographic polynomials approximate their ideal
#' projections over off-axis angles from 0 to 90 deg. Changing normalization
#' does not extend the interval over which those approximations were fitted.
#' The calculation limit `max_off_axis_angle` remains a separate choice.
#'
#' ### Radius in pixels
#'
#' `normalization_radius` must refer to the same angle as the coefficients.
#' To preserve an existing geometry in pixels when changing the reference,
#' multiply its previous normalization radius by \eqn{s}:
#' \deqn{r_{norm,new} = r_{norm,old} \cdot s.}
#'
#' @param projection_coef numeric vector. Finite coefficients of a radial
#'   polynomial evaluated in radians, in increasing order of degree and
#'   without an intercept, as in [set_internal_geometry()]. Must be nonempty.
#' @param normalization_angle numeric vector of length one. Positive off-axis
#'   angle in degrees, at most 180. The input polynomial must have a finite,
#'   positive value at this angle.
#'
#' @return numeric vector of the same length and with the same names as
#'   `projection_coef`, normalized at `normalization_angle`.
#'
#' @seealso [set_internal_geometry()] for the radial model and its validation.
#'
#' @export
#' @examples
#' # Coefficients for each standard projection, normalized at 60 deg.
#' normalize_projection_coef(equidistant_coef, 60)
#' normalize_projection_coef(equal_area_coef, 60)
#' normalize_projection_coef(orthographic_coef, 60)
#'
#' # Preserve a geometry whose radius at 90 deg was 50 pixels.
#' angle <- 60
#' scale <- calc_relative_radius(angle, equal_area_coef)
#' geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(50, 50),
#'   parameters = normalize_projection_coef(equal_area_coef, angle),
#'   normalization_radius = 50 * scale,
#'   normalization_angle = angle
#' )
#' g <- set_internal_geometry(c(100, 100), geometry)
normalize_projection_coef <- function(projection_coef, normalization_angle) {
  .check_vector(projection_coef, "numeric")
  if (!length(projection_coef)) {
    stop("`projection_coef` cannot be empty.", call. = FALSE)
  }
  .check_vector(normalization_angle, "numeric", 1, sign = "positive")
  if (normalization_angle > 180) {
    stop("`normalization_angle` must not exceed 180 degrees.", call. = FALSE)
  }
  reference <- .calc_relative_radius(normalization_angle * pi / 180,
                                      projection_coef)
  if (!is.finite(reference) || reference <= 0) {
    stop("The radial polynomial must have a finite positive value at ",
         "`normalization_angle`.", call. = FALSE)
  }
  coefficients <- projection_coef / reference
  if (any(!is.finite(coefficients))) {
    stop("Normalization produced non-finite coefficients.", call. = FALSE)
  }
  coefficients
}
