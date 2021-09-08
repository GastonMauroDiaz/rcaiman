#' Relative radius image
#'
#' Built a single layer image with relative radius values.
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
#' @param angle Numeric vector. Angles in degrees.
#' @param lens_coef Numeric vector. Polynomial coefficients
#'   of the lens projection function.
#'
#' @noRd
calc_relative_radius <- function(angle, lens_coef) {

  angle <- degree2radian(angle)

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
  round(unname(relative_radius), 2)
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
#' @examples
#' zenith_image(1490, lens("Nikon_FCE9"))
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


#' Azimuth image
#'
#' Built a single layer image with azimuth angles values.
#'
#' @param z \code{\linkS4class{RasterLayer}} built with
#'   \code{\link{zenith_image}}.
#'
#'
#' @return \code{\linkS4class{RasterLayer}}.
#' @export
#'
#' @examples
#' z <- zenith_image(1490, lens("Nikon_FCE9"))
#' azimuth_image(z)
azimuth_image <- function (z)
{
  mask <- is.na(z)

  xy <- xyFromCell(z, seq(length = ncell(z)))
  v <- values(z)
  sph <- geometry::cart2sph(
    matrix(c(xy[, 1] - ncol(z) / 2, xy[, 2] - ncol(z) / 2, values(z)), ncol = 3)
  )

  values(z) <- sph[, 1] * 180 / pi
  values(z) <- values(abs(t(z) - 180)) # to orient North up and West left

  z[mask] <- NA

  z
}
