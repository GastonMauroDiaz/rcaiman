#' Relative radius image
#'
#' Build a single layer image with relative radius values.
#'
#' @param diameter A numeric vector of length one. Diameter in pixels.
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

  # dis2 <- terra::flip(dis1)
  # dis2 <- terra::extend(dis2, terra::ext(0, diameter/2, 0, diameter))
  # terra::ext(dis1) <- terra::ext(0, diameter/2, diameter/2, diameter)
  # dis1 <- terra::resample(dis1, dis2)
  # dis1 <- terra::cover(dis1, dis2)
  # dis2 <- terra::rev(dis1)
  # dis1 <- terra::extend(dis1, terra::ext(diameter/2,
  #                                        diameter,
  #                                        0,
  #                                        diameter))
  # terra::ext(dis2) <- terra::ext(diameter/2, diameter, 0, diameter)
  # dis2 <- terra::resample(dis2, dis1)
  # terra::cover(dis1, dis2)

  # r <- terra::rast(ncol = diameter, nrow = diameter)
  # terra::crs(r) <- "epsg:7589" # https://spatialreference.org/ref/sr-org/7589/
  # terra::ext(r) <- terra::ext(0, diameter, 0, diameter)
  # zenith <- diameter / 2
  # i <- terra::cellFromRowCol(r, zenith, zenith)
  # p1 <- terra::vect(matrix(c(zenith, zenith), ncol = 2), crs = terra::crs(r))
  # dis <- terra::distance(r, p1)
  # i <- terra::cellFromRowCol(r, zenith, diameter)
  # dis <- dis / as.numeric(dis[i])
  # dis[dis > 1] <- NA
  # dis
}

#' Calculate relative radius
#'
#' Calculate the relative radius given a zenith angle and lens function.
#'
#' @param angle Numeric vector. Zenith angles in degrees.
#' @param lens_coef Numeric vector. Polynomial coefficients
#'   of the lens projection function.
#'
#' @noRd
calc_relative_radius <- function(angle, lens_coef) {

  angle <- .degree2radian(angle)

  temp <- cbind(lens_coef, seq(1, length(lens_coef)))
  for (i in 1:length(lens_coef)) {
    if (i == 1) {
      ma <- temp[i, 1] * angle^temp[i, 2]
    } else {
      ma <- rbind(ma, temp[i, 1] * angle^temp[i, 2])
    }
  }

  if (length(lens_coef) == 1) {
    relative_radius <- ma
  } else {
    relative_radius <- apply(ma, 2, sum)
  }
  unname(relative_radius)
}

#' Zenith image
#'
#' Built a single layer image with zenith angles values.
#'
#' @param diameter Numeric vector of length one. Diameter in pixels.
#' @param lens_coef Numeric vector. Polynomial coefficients of the lens
#'   projection function.
#'
#' @return An object of class \linkS4class{SpatRaster} of zenith angles in
#'   degrees, showing a complete hemispherical view, with the zenith on the
#'   center.
#' @export
#'
#' @family Lens functions
#'
#' @examples
#' z <- zenith_image(1490, lens("Nikon_FCE9"))
#' plot(z)
zenith_image <- function (diameter, lens_coef)
{
  # Assign zenith angle by inverting relative radius(R)
  # with a Look Up Table (LUT).
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
  z
}
