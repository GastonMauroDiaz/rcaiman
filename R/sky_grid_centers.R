#' Generate the centers of a sky grid
#'
#' @description
#' Generate external viewing directions at the angular midpoints of the
#' zenith rings and geographic-azimuth sectors used by
#' [sky_grid_segmentation()].
#'
#' @details
#' The grid covers zenith angles from 0 to 90 deg and a full turn in azimuth,
#' measured from north toward east. The reference is the sky, independently
#' of the optical axis or acquisition orientation.
#'
#' Each direction is the midpoint of its two angular intervals, rather than
#' an area centroid. Equal angular widths do not give equal solid angles:
#' the cells occupy less sky near zenith than near the horizon.
#'
#' The centers correspond to the default `first_ring_different = FALSE` in
#' [sky_grid_segmentation()]. Generation is independent of image coverage,
#' raster dimensions, and projection model. Use [reproject_points()] with
#' the destination geometry and orientation to locate these directions in an
#' image, then [indices_from_coord()] or [add_indices()] to assign them to cells.
#'
#' @param angle_width positive numeric vector of length one. Width of both
#'   angular intervals, in degrees. Must divide 90 into a whole number of
#'   rings (and consequently 360 into a whole number of sectors).
#'
#' @return data frame with `zenith` and `azimuth` in radians. Rows follow
#'   increasing zenith intervals; azimuth increases within each interval.
#' @seealso [sky_grid_segmentation()], [reproject_points()], [add_directions()]
#' @export
#' @examples
#' query_directions <- sky_grid_centers(15)
#' head(query_directions)
#'
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#'  )
#' query_coord <- reproject_points(query_directions, geometry, elevation = 32.5)
sky_grid_centers <- function(angle_width) {
  n <- .sky_interval_count(angle_width, 90)
  zenith <- (seq_len(n) - .5) * angle_width * pi / 180
  azimuth <- (seq_len(4 * n) - .5) * angle_width * pi / 180
  directions <- expand.grid(azimuth = azimuth, zenith = zenith,
                              KEEP.OUT.ATTRS = FALSE)
  directions[c("zenith", "azimuth")]
}
