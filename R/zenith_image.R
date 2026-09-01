#' Relative radius image
#'
#' Build a single-layer image with relative radius values
#'
#' @inheritParams zenith_image
#'
#' @noRd
relative_radius_image <- function (diameter, na_below_horizon = TRUE){
  r <- terra::rast(ncol = diameter/2, nrow = diameter/2)
  r <- set_rcaiman_geometry(r)
  zenith <- diameter / 2
  p1 <- terra::vect(matrix(c(zenith, zenith), ncol = 2), crs = terra::crs(r))
  dis <- terra::distance(r, p1)
  dis <- dis / as.numeric(dis[1])
  if (na_below_horizon) dis[dis > 1] <- NA
  dis
}

.quarter_circle_to_full <- function(third_cuadrant) {
  diameter <- ncol(third_cuadrant) * 2
  z1 <- terra::flip(third_cuadrant)
  z2 <- terra::rev(third_cuadrant)
  z4 <- terra::rev(z1)
  terra::ext(z1) <- terra::ext(0, diameter/2, diameter/2, diameter)
  terra::ext(z2) <- terra::ext(diameter/2, diameter, diameter/2, diameter)
  terra::ext(z4) <- terra::ext(diameter/2, diameter, 0, diameter/2)
  z1 <- terra::extend(z1, terra::ext(0, diameter, 0, diameter))
  z2 <- terra::extend(z2, terra::ext(0, diameter, 0, diameter))
  third_cuadrant <- terra::extend(third_cuadrant, terra::ext(0, diameter, 0, diameter))
  z4 <- terra::extend(z4, terra::ext(0, diameter, 0, diameter))
  sum(z1, z2, third_cuadrant, z4, na.rm = TRUE)
}

.zenith_image <- function(diameter, projection_coef, na_below_horizon) {
  x <- relative_radius_image(diameter, na_below_horizon)
  if (na_below_horizon) {
    angle <- seq(0, 90, length.out = nrow(x) + 1)
    R <- calc_relative_radius(angle, projection_coef, units = "degree")
    rcl <- matrix(c(c(0, R[-length(R)]), R, angle), ncol = 3)
  } else {
    angle <- seq(0, 180, length.out = 2 * nrow(x) + 2)
    R <- calc_relative_radius(angle, projection_coef, units = "degree")
    rcl <- matrix(c(c(0, R[-length(R)]), R, angle), ncol = 3)
  }
  z3 <- terra::classify(x, rcl)
  .quarter_circle_to_full(z3)
}

#' Build zenith image
#'
#' @description
#' Build a single-layer image with zenith angle values.
#'
#' @details
#' This implementation assumes a radially symmetric projection model
#' where the zenith corresponds to the center of a circular image
#' and the horizon is represented by a circle at 90° zenith angle.
#'
#' The relative radius is the radial distance expressed as a fraction of the
#' maximum radius (i.e., the horizon). The mapping between zenith angle and
#' relative radius is defined by a polynomial model evaluated in radians (see
#' [calc_relative_radius()]).
#'
#' The equidistant projection model is the standard reference (see
#' [equidistant_coef]). Real lenses deviate from ideal projections. Following
#' [Hemisfer software](https://www.schleppi.ch/patrick/hemisfer/), deviations
#' are modeled using polynomial functions fitted to angular values in radians.
#'
#' @param diameter numeric vector of length one. Diameter in pixels
#'   expressed as an even integer.
#'
#' @inheritParams calc_relative_radius
#'
#' @return [terra::SpatRaster-class] with zenith angles in degrees,
#'   representing a complete hemispherical view with the zenith at the
#'   center. The object carries attribute `projection_coef`.
#'
#' @export
#'
#' @examples
#' z <- zenith_image(600,equidistant_coef)
#' plot(z)
zenith_image <- function (diameter, projection_coef)
{
  .check_vector(diameter, "even_integerish", 1, sign = "positive")
  .check_vector(projection_coef, "numeric", sign = "any")
  validate_projection_function(projection_coef, tol = 1e-2)

  z <- .zenith_image(diameter, projection_coef, na_below_horizon = TRUE)
  names(z) <- "Zenith image"
  attr(z, "projection_coef") <- projection_coef
  .assert_zenith_image(z)
}
