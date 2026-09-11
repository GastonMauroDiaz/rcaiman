#' Reproject points to continuous image coordinates
#'
#' @description
#' Recompute continuous image coordinates after changing the radial projection
#' model or orientation that defines a target image geometry.
#'
#' @details
#' \pkg{rcaiman} represents hemispherical images using a radial projection
#' system, where each point is defined by its angular coordinates (zenith and
#' azimuth).
#'
#' The radial projection model defines how zenith angles are transformed into
#' normalized radial distances (see [calc_relative_radius()]). These distances
#' are then scaled to the continuous image coordinates of a square target with
#' the requested `diameter`.
#'
#' Reprojection preserves angular coordinates. By default, the target geometry
#' has North at the top. A nonzero `orientation` can restore the orientation of
#' a geometry that was previously standardized.
#'
#' Existing continuous image coordinates are replaced. Derived columns tied to
#' the source image frame are removed with [keep_coord()]. Angular
#' directions and other columns are preserved.
#'
#' @param sampling_points `data.frame` with columns `"zen"` and `"azi"` in
#'   radians. Existing `"x"` and `"y_down"` columns are replaced by the target
#'   coordinates. Columns `"col"`, `"row"`, `"y"`, `"relative_radius"`, and
#'   `"image_angle"` are removed. All other columns are preserved.
#' @param orientation numeric vector of length one. Orientation of the target
#'   image (deg), using the convention of [azimuth_image()]. The default is zero
#'   for a North-up target.
#'
#' @inheritParams calc_relative_radius
#' @inheritParams zenith_image
#'
#' @return `data.frame` with columns `"x"` and `"y_down"` followed by preserved
#'   angular directions and other frame-independent columns from
#'   `sampling_points`. Coordinates refer to the target image and use the
#'   convention described in [y_down_from_y()].
#'
#' @export
#' @examples
#' z <- zenith_image(1000, equal_area_coef)
#' a <- azimuth_image(z)
#' sampling_indices <- fibonacci_indices(z, a, 10)
#' reprojected_points <- reproject_points(
#'   sampling_indices,
#'   equidistant_coef,
#'   ncol(z)
#' )
#'
#' plot(!is.na(z))
#' points(xy_from_indices(
#'          sampling_indices$col,
#'          sampling_indices$row,
#'          nrow(z)
#'        ),
#'        col = "blue", pch = 10)
#' points(reprojected_points$x,
#'        y_from_y_down(reprojected_points$y_down, nrow(z)),
#'        col = "darkgreen", pch = 10)
reproject_points <- function(sampling_points,
                             projection_coef,
                             diameter,
                             orientation = 0) {

  .check_columns(sampling_points,
                 required_columns = c("zen", "azi"))
  .check_vector(diameter, "even_integerish", 1, sign = "positive")
  .check_vector(orientation, "numeric", 1)

  zen <- sampling_points$zen
  azi <- sampling_points$azi

  R <- .calc_relative_radius(zen, projection_coef)

  orientation <- orientation * pi / 180
  theta <- azi + orientation + pi/2
  r <- R * (pi/2)

  x <- r * cos(theta)
  y <- r * sin(theta)

  n <- diameter
  res <- pi / n

  x <- (x + pi/2) / res
  y_down <- (pi/2 - y) / res

  sampling_points$x <- x
  sampling_points$y_down <- y_down
  keep_coord(sampling_points)
}
