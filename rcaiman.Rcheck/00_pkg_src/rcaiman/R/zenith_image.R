#' Relative radius image
#'
#' Build a single layer image with relative radius values.
#'
#' @param diameter A numeric vector of length one. Diameter in pixels.
#'
#' @noRd
relative_radius_image <- function (diameter)
{
  r <- raster(ncol = diameter, nrow = diameter)
  extent(r) <- c(0, diameter, 0, diameter)
  values(r) <- 1
  projection(r) <- NA

  zenith <- diameter / 2
  dis <- distanceFromPoints(r, matrix(c(zenith, zenith), ncol = 2))
  dis[dis > zenith] <- NA
  dis <- dis / zenith

  values(r) <- values(dis)

  r
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
#'
#' @param diameter Numeric vector of length one. Diameter in pixels.
#' @param lens_coef Numeric vector. Polynomial coefficients
#'   of the lens projection function.
#'
#' @return \code{\linkS4class{RasterLayer}}.
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

  x <- relative_radius_image(diameter)

  angle <- seq(0, 90,  length.out = nrow(x) + 1)
  R <- calc_relative_radius(angle, lens_coef)

  rcl <- matrix(c(c(0, R[-length(R)]), R, angle), ncol = 3)

  reclassify(x, rcl)
}
