#' CIE sky model
#'
#' Written by Gaston Mauro Diaz based on Pascal code by Mait Lang.
#'
#' Angles should be provided in radians.
#'
#' @param AzP Numeric vector. Azimuth angle of a sky point.
#' @param Zp Numeric vector. Zenith Angle of a sky point.
#' @param AzS Numeric vector of length one. Azimuth angle of the sun disc.
#' @param Zs Numeric vector of length one. Zenith angle of the sun disc.
#' @param .a,.b,.c,.d,.e Numeric vector of length one. Sky model parameter.
#'
#' @noRd
#' @references http://dx.doi.org/10.1016/j.energy.2016.02.054
#'
#' @return Numeric vector of length equal to AzP length.
.cie_sky_model <- function(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e) {
  # calculate distances
  ## between sky point and Sun
  chi_s <- calc_spherical_distance(Zp, AzP, Zs, AzS)
  ## between sky point and zenith
  chi_z <- Zs

  # Gradation function
  .Phi <- function(z) 1 + .a * exp(.b / cos(z))

  # Indicatrix function
  .f <- function(chi) 1 + .c * (exp(.d * chi) - exp(.d * pi/2)) + .e * cos(chi)^2

  unname((.Phi(Zp) * .f(chi_s)) / (.Phi(0) * .f(chi_z)))
}


#' CIE sky model raster
#'
#' @inheritParams sky_grid_segmentation
#' @param sun_zenith_azimuth Numeric vector of length two. The solar disk
#'   center represented with zenith and azimuth angles in degrees.
#' @param sky_coef Numeric vector of length five. Parameters of the sky model.
#'
#' @export
#'
#' @examples
#' z <- zenith_image(50, lens())
#' a <- azimuth_image(z)
#' path <- system.file("external", package = "rcaiman")
#' skies <- read.csv(file.path(path, "15_CIE_standard_skies.csv"))
#' # parameters are from http://dx.doi.org/10.1016/j.energy.2016.02.054
#' sky_coef <- skies[4,1:5]
#' sun_zenith_azimuth <- c(45, 0)
#' plot(cie_sky_image(z, a, sun_zenith_azimuth, sky_coef))
cie_sky_image <- function(z, a, sun_zenith_azimuth, sky_coef) {
  .is_single_layer_raster(z)
  .is_single_layer_raster(a)
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  stopifnot(length(sun_zenith_azimuth) == 2)
  stopifnot(length(sky_coef) == 5)

  Zp <- .degree2radian(z[])
  AzP <- .degree2radian(a[])

  Zs <- .degree2radian(sun_zenith_azimuth[1])
  AzS <- .degree2radian(sun_zenith_azimuth[2])

  relative_luminance <- .cie_sky_model(AzP, Zp, AzS, Zs,
                                       as.numeric(sky_coef[1]),
                                       as.numeric(sky_coef[2]),
                                       as.numeric(sky_coef[3]),
                                       as.numeric(sky_coef[4]),
                                       as.numeric(sky_coef[5]))
  terra::values(z) <- relative_luminance
  z[is.infinite(z)] <- 0
  names(z) <- "Relative radiance or luminance"
  z
}
