#' Build a relative solid-angle weight image
#'
#' @description
#' Calculate relative pixel solid-angle weights analytically from a radial
#' projection model and the internal off-axis angles.
#'
#' @param g [terra::SpatRaster-class] from
#'   [set_internal_geometry()], with `off_axis` and `internal_azimuth` layers
#'   in degrees. Only `off_axis` is used; missing azimuths do not affect weights.
#' @param geometry named list following [set_internal_geometry()], with
#'   `model = "radial_projection"`. Supply the geometry used to construct
#'   `g`, including its normalization reference. Anisotropic
#'   models are not supported by this analytical calculation.
#'
#' @details
#' Inclining the optical system changes the external direction of each pixel
#' but not the solid angle it represents. Weights therefore use off-axis
#' angle, not external zenith angle, and can be reused with different
#' acquisition orientations.
#'
#' For the radial polynomial \eqn{f} defined in [set_internal_geometry()],
#' evaluated with off-axis angle \eqn{\theta'} in radians, the weight is
#' \deqn{W(\theta') = \frac{\sin(\theta')}{f(\theta') f'(\theta')}.}
#' At the optical axis the limit is \eqn{W(0) = 1/a_1^2}, where \eqn{a_1}
#' is the first polynomial coefficient. The model must have a positive
#' derivative at the evaluated directions to give finite weights.
#'
#' These are relative weights. Divide by `geometry$normalization_radius^2`
#' to obtain the local solid angle per unit pixel area in steradians. This
#' is a differential approximation evaluated at pixel centers, not an
#' integration over each pixel boundary. The normalization angle need not
#' be 90 deg. Changing that reference rescales the relative weights; dividing
#' by the corresponding squared radius gives the same physical result.
#'
#' Cropping or masking the internal geometry preserves the weights at the
#' remaining pixels. Resampling to a different image scale requires the
#' geometry of that new representation.
#'
#' @return Single-layer [terra::SpatRaster-class] on the same grid as
#'   `g`, with relative solid-angle weights and `NA` where
#'   off-axis angle is missing. No calibration attributes are attached.
#' @seealso [thr_isodata()], [apply_by_direction()], [merge_segments()],
#'   [normalize_projection_coef()]
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
#' w <- area_weight_image(g, geometry)
#' plot(w)
#' solid_angle <- w / geometry$normalization_radius^2
area_weight_image <- function(g, geometry) {
  .check_internal_geometry(geometry)
  if (geometry$model != "radial_projection") {
    stop("`geometry` must use `radial_projection` for analytical weights.",
         call. = FALSE)
  }
  .assert_spatraster(g)
  if (!identical(names(g), c("off_axis", "internal_azimuth"))) {
    stop("`g` must have layers off_axis and internal_azimuth, ",
         "in that order.", call. = FALSE)
  }
  angle <- terra::values(g$off_axis, mat = FALSE)
  .check_vector(angle, "numeric", allow_na = TRUE, name = "off_axis")
  limit <- if ("max_off_axis_angle" %in% names(geometry))
    geometry$max_off_axis_angle else 90
  if (any(angle < 0 | angle > limit, na.rm = TRUE)) {
    stop("Off-axis angles must be within the radial model's calculation limit.",
         call. = FALSE)
  }
  theta <- angle * pi / 180
  coef <- geometry$parameters
  radius <- .calc_relative_radius(theta, coef)
  slope <- rep(0, length(theta))
  for (i in rev(seq_along(coef))) slope <- slope * theta + i * coef[i]
  if (any(slope <= 0, na.rm = TRUE)) {
    stop("The radial model must have a positive derivative at evaluated angles.",
         call. = FALSE)
  }
  weight <- sinpi(angle / 180) / (radius * slope)
  axial <- !is.na(angle) & angle == 0
  weight[axial] <- 1 / coef[1]^2
  if (any(!is.finite(weight[!is.na(angle)]))) {
    stop("The radial model produced non-finite solid-angle weights.", call. = FALSE)
  }
  out <- terra::rast(g$off_axis)
  terra::values(out) <- weight
  names(out) <- "Area weight image"
  out
}
