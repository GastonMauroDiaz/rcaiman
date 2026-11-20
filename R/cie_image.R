#' CIE sky model
#'
#' Written by Gaston Mauro Diaz based on Pascal code by Mait Lang.
#'
#' Angles should be provided in radians.
#'
#' @param AzP numeric vector. Azimuth angle of a sky point.
#' @param Zp numeric vector. Zenith Angle of a sky point.
#' @param AzS numeric vector of length one. Azimuth angle of the sun disc.
#' @param Zs numeric vector of length one. Zenith angle of the sun disc.
#' @param .a,.b,.c,.d,.e numeric vector of length one. Sky model parameter.
#'
#' @noRd
#' @references http://dx.doi.org/10.1016/j.energy.2016.02.054
#'
#' @return numeric vector of length equal to AzP length.
.cie_sky_model <- function(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e) {
  # calculate angular distance between sky point and Sun
  Chi <- calc_spherical_distance(Zp, AzP, Zs, AzS)

  # Gradation function
  Phi_Z <- 1 + .a * exp(.b / cos(Zp))
  Phi_0 <- 1 + .a * exp(.b)
  gradation <- Phi_Z / Phi_0

  # Indicatrix function
  F_Chi <-  1 + .c * (exp(.d * Chi) - exp(.d * pi/2)) + .e * cos(Chi)^2
  F_Zs  <-  1 + .c * (exp(.d *  Zs) - exp(.d * pi/2)) + .e * cos(Zs)^2
  indicatrix <- F_Chi / F_Zs

  unname(gradation * indicatrix)
}
.cie_sky_model <- compiler::cmpfun(.cie_sky_model)



#' CIE sky image
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Evaluate the CIE General Sky model at the external viewing direction of
#' each pixel, producing relative radiance or luminance.
#'
#' @details
#' Coefficient sets and the formulation are available in [cie_table]. The
#' angles in `d` are converted from degrees to radians before evaluation.
#' Directions are external: the acquisition orientation is already represented
#' in `d`, independently of the lens projection and the image grid.
#'
#' The model is normalized by its value at zenith, not at the optical axis or
#' the brightest pixel. This reference is defined by the model even if the
#' image does not include zenith. For a model with a defined, nonzero zenith
#' reference, its relative value at zenith is one. Multiply the output by a
#' separately obtained zenith radiance, luminance, or radiometrically corrected
#' DN to express predictions on that scale. The function does not estimate
#' this reference from the image.
#'
#' Only sky directions with zenith angles from 0 to 90 deg are evaluated.
#' At zenith, azimuth is unnecessary and may be `NA`; elsewhere both angles
#' are required. `m` can further restrict evaluation without changing the
#' model normalization.
#'
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, in that order.
#' @param sun_angles data frame of one row with finite numeric columns
#'   `zenith` and `azimuth` in radians, as returned by [add_directions()] for
#'   external directions. Column order is arbitrary and additional columns
#'   are ignored. The solar direction
#'   must use the same external reference as `d`. Zenith must be nonnegative;
#'   values greater than `96 * pi / 180` produce a warning. Azimuth must be
#'   in `[0, 2 * pi)`.
#' @param sky_coef finite numeric vector of length five. CIE model coefficients
#'   in the order `a`, `b`, `c`, `d`, `e`, from [cie_table] or [fit_cie_model()].
#' @param m optional single-layer logical, NA-free [terra::SpatRaster-class]
#'   on the same grid as `d`. Only cells marked `TRUE` are evaluated. With
#'   `NULL`, all available sky directions are evaluated.
#'
#' @return Single-layer numeric [terra::SpatRaster-class] on the grid of `d`,
#'   named `Relative radiance or luminance`. Values are relative luminance for
#'   the standard sky coefficients, or relative radiance or luminance according
#'   to the data used to fit `sky_coef`. Directions below the horizon, outside
#'   the mask, or missing a required angle receive `NA`.
#'
#' @note
#' Negative and nonfinite model predictions are returned as `NA`. Supplying
#' custom coefficients does not trigger fitting or coefficient adjustment.
#'
#' @seealso [cie_table], [fit_cie_model()], [direction_image()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' sky_coef <- as.numeric(cie_table[4, 1:5])
#' # Cell containing the marked solar disc; extracted angles are in radians.
#' sun_angles <- add_directions(data.frame(col = 40, row = 30), d = d)
#' # Alternatively: data.frame(zenith = 45 * pi/180, azimuth = 10 * pi/180)
#' sky <- cie_image(d, sun_angles, sky_coef)
#' plot(sky)
#'
#' m <- !is.na(d$zenith) & d$zenith < 80
#' sky <- cie_image(d, sun_angles, sky_coef, m = m)
#' plot(sky)
#'
#' # A partial field of view need not include the normalization direction.
#' side_geometry <- geometry
#' side_geometry$normalization_radius <- 100
#' side_g <- set_internal_geometry(c(80, 60), side_geometry)
#' side_d <- direction_image(side_g, elevation = 32.5)
#' side_sky <- cie_image(side_d, sun_angles, sky_coef)
#' plot(side_sky)
cie_image <- function(d, sun_angles, sky_coef, m = NULL) {
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  .check_cie_sun_angles(sun_angles)
  .check_vector(sky_coef, "numeric", 5, sign = "any")

  zenith <- terra::values(z, mat = FALSE)
  azimuth <- terra::values(a, mat = FALSE)
  at_zenith <- !is.na(zenith) & zenith == 0
  azimuth[at_zenith] <- 0
  selected <- !is.na(zenith) & zenith <= 90 & !is.na(azimuth)
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
    selected <- selected & terra::values(m, mat = FALSE)
  }

  relative <- rep(NA_real_, terra::ncell(z))
  if (any(selected)) {
    predicted <- .cie_sky_model(
      .degree2radian(azimuth[selected]), .degree2radian(zenith[selected]),
      sun_angles[["azimuth"]], sun_angles[["zenith"]],
      sky_coef[1], sky_coef[2], sky_coef[3], sky_coef[4], sky_coef[5]
    )
    predicted[!is.finite(predicted) | predicted < 0] <- NA_real_
    relative[selected] <- predicted
  }
  sky_cie <- terra::rast(z)
  terra::values(sky_cie) <- relative
  names(sky_cie) <- "Relative radiance or luminance"
  sky_cie
}

.check_cie_sun_angles <- function(sun_angles) {
  if (!is.data.frame(sun_angles) || nrow(sun_angles) != 1L ||
      !all(c("zenith", "azimuth") %in% names(sun_angles))) {
    stop("`sun_angles` must be a data frame of one row with numeric columns ",
         "'zenith' and 'azimuth' in radians.", call. = FALSE)
  }
  .check_vector(sun_angles[["zenith"]], "numeric", 1, sign = "nonnegative")
  .check_vector(sun_angles[["azimuth"]], "numeric", 1, sign = "nonnegative")
  if (sun_angles[["azimuth"]] >= 2 * pi) {
    stop("Sun disk azimuth angle must lie within [0, 2 * pi) radians.",
         call. = FALSE)
  }
  if (sun_angles[["zenith"]] > 96 * pi / 180) {
    warning("Sun disk zenith angle should not be greater than 96 deg ",
            "(96 * pi / 180 radians).", call. = FALSE)
  }
  invisible(TRUE)
}
