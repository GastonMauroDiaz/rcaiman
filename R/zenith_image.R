#' Relative radius image
#'
#' Build a single-layer image with relative radius values
#'
#' @inheritParams zenith_image
#'
#' @noRd
relative_radius_image <- function (diameter)
{
  r <- terra::rast(ncol = diameter/2, nrow = diameter/2)
  terra::crs(r) <- "epsg:7589" # https://spatialreference.org/ref/sr-org/7589/
  terra::ext(r) <- terra::ext(0, diameter/2, 0, diameter/2)
  zenith <- diameter / 2
  p1 <- terra::vect(matrix(c(zenith, zenith), ncol = 2), crs = terra::crs(r))
  dis <- terra::distance(r, p1)
  dis <- dis / as.numeric(dis[1])
  dis[dis > 1] <- NA
  dis
}



#' Build Zenith image
#'
#' Build a single layer-image with zenith angle values, assuming upwards-looking
#' hemispherical photography with the optical axis vertically aligned.
#'
#'
#' @param diameter Numeric vector of length one. Diameter in pixels expressed as
#'   an even integer. The latter is to simplify calculations by having the
#'   zenith point located between pixels. Snapping the zenith point between
#'   pixels does not affect accuracy because half-pixel is less than the
#'   uncertainty in localizing the circle within the picture.
#' @param lens_coef Numeric vector. Polynomial coefficients of the lens
#'   projection function. See [calibrate_lens()].
#'
#' @return An object of class [SpatRaster-class] of zenith angles in degrees,
#'   showing a complete hemispherical view with the zenith on the center.
#' @export
#'
#' @examples
#' z <- zenith_image(600, lens("Nikon_FCE9"))
#' plot(z)
zenith_image <- function (diameter, lens_coef)
{
  # Assign zenith angle by inverting relative radius(R) with a Look Up Table
  stopifnot(.is_whole(diameter))
  stopifnot(.is_even(diameter))

  x <- relative_radius_image(diameter)
  angle <- seq(0, 90,  length.out = nrow(x) + 1)
  R <- calc_relative_radius(angle, lens_coef)
  rcl <- matrix(c(c(0, R[-length(R)]), R, angle), ncol = 3)
  z3 <- terra::classify(x, rcl)
  z1 <- terra::flip(z3)
  z2 <- terra::rev(z3)
  z4 <- terra::rev(z1)
  terra::ext(z1) <- terra::ext(0, diameter/2, diameter/2, diameter)
  terra::ext(z2) <- terra::ext(diameter/2, diameter, diameter/2, diameter)
  terra::ext(z4) <- terra::ext(diameter/2, diameter, 0, diameter/2)
  z1 <- terra::extend(z1, terra::ext(0, diameter, 0, diameter))
  z2 <- terra::extend(z2, terra::ext(0, diameter, 0, diameter))
  z3 <- terra::extend(z3, terra::ext(0, diameter, 0, diameter))
  z4 <- terra::extend(z4, terra::ext(0, diameter, 0, diameter))
  z <- sum(z1, z2, z3, z4, na.rm = TRUE)
  names(z) <- "Zenith image"
  z
}
