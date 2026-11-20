#' Reproject directions to a radial destination
#'
#' @description
#' Map viewing directions to continuous image coordinates in a radial
#' representation with a specified orientation.
#'
#' @param sampling_points `data.frame` from [add_directions()], with numeric
#'   columns `"zenith"` and `"azimuth"` in radians. Zenith is in `[0, pi]`
#'   and azimuth in `[0, 2 * pi)`. Missing angles may be `NA`.
#' @param geometry named list following [set_internal_geometry()], with
#'   `model = "radial_projection"`. Describes the destination, including its
#'   center of vision, projection coefficients, normalization radius and angle,
#'   and optional `max_off_axis_angle` (90 deg when absent). Additional fields
#'   and list classes are accepted.
#' @param sensor_rotation numeric vector of length one. Clockwise rotation of
#'   the destination sensor in degrees, using the convention of
#'   [direction_image()]. Defaults to zero (landscape).
#' @param azimuth_rotation numeric vector of length one. Horizontal turn of
#'   the destination to the right from south, in degrees, as in
#'   [direction_image()]. Defaults to zero.
#' @param elevation numeric vector of length one. Elevation of the destination
#'   optical axis above the horizontal in degrees, in `[-90, 90]`, as in
#'   [direction_image()]. Defaults to 90 (zenith).
#'
#' @details
#' The orientation arguments describe the destination; source orientation is
#' already expressed by `zenith` and `azimuth`. Their defaults place zenith at
#' the center and North at the top, with East to the left and West to the right.
#'
#' Directions are first expressed in the destination's internal frame with
#' [internal_directions_from_external()], then mapped to image coordinates by
#' [coord_from_internal_directions()]. The radial polynomial
#' converts the resulting off-axis angle to relative radius. Multiplication
#' by `normalization_radius` and the internal azimuth then locate the point
#' around `center_of_vision`.
#'
#' No interpolation or snapping is applied. Existing `x` and `y_down` are
#' replaced; derived columns tied to the previous frame are removed with
#' [keep_coord()]. Angular coordinates, other values, and row order are
#' preserved, including rows without a position in the destination.
#'
#' Directions beyond `max_off_axis_angle`, measured from the destination's
#' optical axis, have `NA` target coordinates. Missing directions also have
#' `NA` coordinates, except that zenith and nadir are defined even when their
#' geographic azimuth is `NA`. The destination's optical-axis direction maps
#' to `center_of_vision`; the opposite direction has no unique radial position and
#' receives `NA` coordinates, even when the angular limit is 180 deg.
#' No image dimensions are needed, so coordinates are not clipped to a frame.
#'
#' @return `data.frame` with `"x"` and `"y_down"` followed by the preserved
#'   columns. Coordinates follow [y_down_from_y()].
#'
#' @seealso [add_directions()], [reproject_image()], [set_internal_geometry()]
#'
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(50, 50),
#'   parameters = equidistant_coef,
#'   normalization_radius = 50,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(100, 100), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' sampling_indices <- data.frame(col = c(40, 50, 60), row = c(45, 50, 55))
#' sampling_indices <- rem_unselected(sampling_indices, !is.na(d[["zenith"]]))
#' points <- add_directions(sampling_indices, d = d)
#' reproject_points(points, geometry)
#' # The same geographic directions in a destination inclined toward west.
#' reproject_points(points, geometry, sensor_rotation = 90,
#'                  azimuth_rotation = 90, elevation = 32.5)
reproject_points <- function(sampling_points, geometry,
                             sensor_rotation = 0, azimuth_rotation = 0,
                             elevation = 90) {
  .check_radial_destination(geometry)
  internal <- internal_directions_from_external(
    sampling_points, sensor_rotation, azimuth_rotation, elevation)
  coord <- coord_from_internal_directions(internal, geometry)
  sampling_points$x <- coord$x
  sampling_points$y_down <- coord$y_down
  keep_coord(sampling_points)
}

.check_radial_destination <- function(geometry) {
  .check_internal_geometry(geometry)
  if (geometry[["model"]] != "radial_projection") {
    stop("The destination `geometry` must use `radial_projection`.",
         call. = FALSE)
  }
  invisible(TRUE)
}
