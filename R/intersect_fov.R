#' Intersect two fields of view
#'
#' @description
#' Compute the common angular region of two field-of-view contours,
#' independently of an image grid.
#'
#' @details
#' Consecutive directions are joined by the shorter great-circle arc, including
#' the last and first vertices. Intersections of these arcs can add vertices to
#' the result. Computation uses spherical polygons through \pkg{s2}; no image
#' reprojection, convex hull, or contour simplification is applied.
#'
#' Vertex order determines the interior on the sphere. Retain the order supplied
#' by [extract_fov()]: clockwise around the covered region in the
#' image, with the usual internal-azimuth convention. Reversing this order selects
#' the complementary spherical region. This convention also supports fields of
#' view extending beyond 90 deg off-axis.
#'
#' Both contours must describe directions in the same internal reference frame.
#' For cross-calibration, acquisition must align the optical axes and orient the
#' tops of the sensors identically. The function does not estimate this alignment.
#'
#' This operation supports a single, simple contour with nonzero spherical area.
#' Self-intersections and retraced edges do not define an admissible contour.
#' An empty table represents an empty field of view. An intersection requiring
#' multiple contours, including holes, produces an error. Contact along an edge
#' or at a point alone has no area and returns an empty table.
#'
#' The contour extracted from an image follows pixel centers, as explained in
#' [extract_fov()]. Intersecting contours does not recover coverage
#' beyond that boundary. Use [clip_fov()] to represent the result on
#' an image grid, or [select_directions_in_fov()] to select angular queries.
#'
#' @param fov_1,fov_2 data frames with numeric columns
#'   `off_axis` in `[0, pi]` and `internal_azimuth` in `[0, 2 * pi)`, in radians,
#'   as returned by [extract_fov()]. Vertices must be ordered along
#'   the contour. The closing vertex may be repeated. Missing angles are not
#'   allowed. Additional columns are ignored.
#'
#' @return data frame with `off_axis` and `internal_azimuth` in radians. The
#'   ordered contour repeats its first direction at the end. Its starting
#'   vertex need not match either input. An empty intersection has zero rows.
#'
#' @seealso [extract_fov()], [select_directions_in_fov()],
#'   [clip_fov()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(20, 15),
#'   normalization_radius = 25,
#'   normalization_angle = 90
#' )
#' uncal_g <- set_internal_geometry(c(40, 30), geometry)
#' geometry$center_of_vision <- c(15, 20)
#' cal_g <- set_internal_geometry(c(30, 40), geometry)
#' uncal_fov <- extract_fov(uncal_g)
#' cal_fov <- extract_fov(cal_g)
#' common_fov <- intersect_fov(uncal_fov, cal_fov)
#' common_uncal_g <- clip_fov(uncal_g, common_fov)
#' common_cal_g <- clip_fov(cal_g, common_fov)
#' query_directions <- fibonacci_directions(10, max_off_axis_angle = 80)
#' query_directions <- select_directions_in_fov(query_directions, common_fov)
intersect_fov <- function(fov_1, fov_2) {
  first <- .spherical_field_of_view(fov_1, "fov_1")
  second <- .spherical_field_of_view(fov_2, "fov_2")
  common <- s2::s2_intersection(first, second,
    options = s2::s2_options(model = "closed", dimensions = "polygon"))
  .field_of_view_from_sphere(common)
}
