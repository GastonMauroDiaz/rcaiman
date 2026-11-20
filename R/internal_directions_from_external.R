#' Convert between external and internal viewing directions
#'
#' @description
#' Express viewing directions in the instrument's reference frame or in the
#' external frame, using the declared acquisition orientation.
#'
#' @details
#' These functions express the direction of a viewing ray in either the
#' instrument's frame or the external angular frame defined by the local
#' vertical and north. Both descriptions refer to the same viewing ray.
#'
#' They do not need a projection model, image dimensions, or raster.
#' All input and output direction
#' columns are in radians; the three orientation arguments are in degrees, as
#' in [direction_image()]. The same orientation describes both conversions.
#'
#' `external_directions_from_internal()` applies the sensor rotation, elevation,
#' and horizontal turn defined by [direction_image()].
#' `internal_directions_from_external()` reverses that transformation. The latter
#' does not ask for negated orientation arguments: supply the actual orientation
#' of the instrument in which the internal directions are to be expressed.
#'
#' Internal directions use `off_axis` and clockwise `internal_azimuth`, measured
#' from the sensor top. External directions use `zenith` and geographic `azimuth`.
#' In the default upward orientation with sensor top toward north, off-axis and
#' zenith angles coincide, while the two azimuths have opposite signs modulo
#' `2*pi`. No hemisphere or field-of-view restriction is applied.
#'
#' At either pole of the input reference frame, the polar angle alone defines
#' the direction and its azimuth may be `NA`. Otherwise, a missing angle produces two
#' missing output angles. At a pole of the output reference frame, the polar
#' angle is zero or pi and the azimuth is `NA`. Numerical residuals within
#' `16 * .Machine$double.eps` in the transverse component of a unit vector are
#' treated as a pole, consistently with [direction_image()].
#'
#' Input columns are retained. Only the two destination-angle columns are added
#' or replaced. Columns describing coordinates, cell indices, or the other
#' reference frame are not recomputed. Each conversion reads only its source
#' angle columns; additional columns do not influence the transformation.
#'
#' To locate an external direction in an image, first convert it to the
#' instrument's internal frame, then use [coord_from_internal_directions()].
#' That second operation currently supports radial projection. The combined
#' operation is also available through [reproject_points()]. Use
#' [internal_directions_from_coord()] for the reverse path from continuous
#' coordinates, or [add_directions()] to extract directions from raster cells.
#'
#' @param directions data frame containing the source angles in radians.
#'   For `internal_directions_from_external()`, supply numeric `zenith` in
#'   `[0, pi]` and `azimuth` in `[0, 2*pi)`. For
#'   `external_directions_from_internal()`, supply numeric `off_axis` in
#'   `[0, pi]` and `internal_azimuth` in `[0, 2*pi)`. Missing angles may be `NA`.
#'   Additional columns are preserved.
#' @inheritParams direction_image
#'
#' @return Data frame preserving input row order, row names, repeated rows, and
#'   other columns. `internal_directions_from_external()` adds or replaces
#'   `off_axis` and `internal_azimuth`. `external_directions_from_internal()`
#'   adds or replaces `zenith` and `azimuth`. Both output angles are in radians.
#'
#' @seealso [direction_image()], [add_directions()],
#'   [coord_from_internal_directions()], [select_directions_in_fov()]
#' @export
#' @examples
#' geometry <- list(model = "radial_projection", parameters = equidistant_coef,
#'   center_of_vision = c(40, 30), normalization_radius = 30,
#'   normalization_angle = 90)
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, azimuth_rotation = 60, elevation = 83)
#'
#' # Locate a known external direction in this instrument's image.
#' sun_angles <- data.frame(zenith = 45 * pi/180, azimuth = 10 * pi/180)
#' internal <- internal_directions_from_external(
#'   sun_angles, azimuth_rotation = 60, elevation = 83)
#' sun_coord <- coord_from_internal_directions(internal, geometry)
#' sun_indices <- add_indices(sun_coord)
#'
#' # Recover external directions directly from continuous image coordinates.
#' internal <- internal_directions_from_coord(sun_coord, geometry)
#' external <- external_directions_from_internal(
#'   internal, azimuth_rotation = 60, elevation = 83)
#' external[c("zenith", "azimuth")]
#'
#' # Select external queries within the instrument's field of view.
#' queries <- sky_grid_centers(15)
#' internal <- internal_directions_from_external(
#'   queries, azimuth_rotation = 60, elevation = 83)
#' selected <- select_directions_in_fov(internal, extract_fov(g))
#' selected[c("zenith", "azimuth")]
internal_directions_from_external <- function(directions,
                                               sensor_rotation = 0,
                                               azimuth_rotation = 0,
                                               elevation = 90) {
  .check_direction_orientation(sensor_rotation, azimuth_rotation, elevation)
  .check_columns(directions, required_columns = c("zenith", "azimuth"))
  .check_vector(directions$zenith, "numeric", allow_na = TRUE, name = "zenith")
  .check_vector(directions$azimuth, "numeric", allow_na = TRUE, name = "azimuth")
  if (any(directions$zenith < 0 | directions$zenith > pi, na.rm = TRUE)) {
    stop("`zenith` must be in [0, pi] radians.", call. = FALSE)
  }
  if (any(directions$azimuth < 0 | directions$azimuth >= 2*pi, na.rm = TRUE)) {
    stop("`azimuth` must be in [0, 2*pi) radians.", call. = FALSE)
  }
  angles <- .internal_direction_angles(as.matrix(directions[c("zenith", "azimuth")]),
    sensor_rotation, azimuth_rotation, elevation)
  directions$off_axis <- angles[, 1]
  directions$internal_azimuth <- angles[, 2]
  directions
}

#' @rdname internal_directions_from_external
#' @export
external_directions_from_internal <- function(directions,
                                               sensor_rotation = 0,
                                               azimuth_rotation = 0,
                                               elevation = 90) {
  .check_direction_orientation(sensor_rotation, azimuth_rotation, elevation)
  .check_internal_directions(directions, "directions", allow_na = TRUE)
  angles <- .external_direction_angles(as.matrix(directions[c("off_axis", "internal_azimuth")]),
    sensor_rotation, azimuth_rotation, elevation)
  directions$zenith <- angles[, 1]
  directions$azimuth <- angles[, 2]
  directions
}

.external_direction_angles <- function(angles, sensor_rotation, azimuth_rotation, elevation) {
  # Rotation of unit rays, starting with clockwise azimuth from sensor top.
  axial <- !is.na(angles[, 1]) & angles[, 1] %in% c(0, pi)
  angles[axial, 2] <- 0
  theta <- angles[, 1] / pi
  alpha <- angles[, 2] / pi + (sensor_rotation %% 360) / 180
  right <- sinpi(theta) * sinpi(alpha)
  top <- sinpi(theta) * cospi(alpha)
  forward <- cospi(theta)
  sin_elevation <- sinpi(elevation / 180)
  cos_elevation <- cospi(elevation / 180)
  horizontal <- forward * cos_elevation - top * sin_elevation
  up <- forward * sin_elevation + top * cos_elevation
  yaw <- (azimuth_rotation %% 360) / 180
  east <- -right * cospi(yaw) - horizontal * sinpi(yaw)
  north <- right * sinpi(yaw) - horizontal * cospi(yaw)
  .angles_from_ray(east, north, up, c("zenith", "azimuth"))
}

.internal_direction_angles <- function(angles, sensor_rotation, azimuth_rotation, elevation) {
  # Geographic east/north/up rays, including poles with unspecified azimuth.
  theta <- angles[, 1] / pi
  alpha <- angles[, 2] / pi
  east <- sinpi(theta) * sinpi(alpha)
  north <- sinpi(theta) * cospi(alpha)
  up <- cospi(theta)
  vertical <- !is.na(angles[, 1]) & angles[, 1] %in% c(0, pi)
  east[vertical] <- north[vertical] <- 0
  # Undo horizontal turn, elevation, and sensor rotation.
  yaw <- (azimuth_rotation %% 360) / 180
  right <- -east * cospi(yaw) + north * sinpi(yaw)
  horizontal <- -east * sinpi(yaw) - north * cospi(yaw)
  forward <- horizontal * cospi(elevation / 180) + up * sinpi(elevation / 180)
  top <- -horizontal * sinpi(elevation / 180) + up * cospi(elevation / 180)
  roll <- (sensor_rotation %% 360) / 180
  internal_right <- right * cospi(roll) - top * sinpi(roll)
  internal_top <- right * sinpi(roll) + top * cospi(roll)
  .angles_from_ray(internal_right, internal_top, forward,
    c("off_axis", "internal_azimuth"))
}

.angles_from_ray <- function(right, top, forward, columns) {
  transverse <- sqrt(right^2 + top^2)
  polar <- atan2(transverse, forward)
  azimuth <- atan2(right, top) %% (2*pi)
  axial <- !is.na(transverse) & transverse <= 16 * .Machine$double.eps
  polar[axial] <- ifelse(forward[axial] >= 0, 0, pi)
  azimuth[axial] <- NA_real_
  azimuth[!is.na(azimuth) & azimuth >= 2*pi] <- 0
  result <- cbind(polar, azimuth)
  colnames(result) <- columns
  result
}
