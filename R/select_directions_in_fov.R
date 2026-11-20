#' Select internal directions within a field of view
#'
#' @description
#' Retain directions belonging to the angular region delimited by a
#' field-of-view contour.
#'
#' @details
#' Membership is computed on the sphere, including the contour boundary.
#' Boundary comparisons allow `64 * .Machine$double.eps` radians for numerical
#' roundoff when converting directions.
#' See [intersect_fov()] for the interpretation of contour order
#' and great-circle edges.
#'
#' Input directions are used directly, without snapping to pixels. Additional
#' columns, including `col` and `row`, do not affect selection. A row with a
#' missing angle is excluded. Angles and other columns in retained rows are
#' not modified.
#'
#' For paired measurements for use in [rcaiman.registry::sample_radiometry()],
#' intersect the fields of view and filter the query table, tipically obtained
#' with [fibonacci_directions()], with the common contour. Supply the resulting
#' table to both sampling operations.
#'
#' @param directions data frame with numeric `off_axis` in `[0, pi]` and
#'   `internal_azimuth` in `[0, 2 * pi)`, both in radians, as returned by
#'   [fibonacci_directions()] or [add_directions()]. Missing angles may be `NA`.
#'   Additional columns are preserved.
#' @param fov data frame representing an ordered contour, as returned
#'   by [extract_fov()] or [intersect_fov()]. Columns
#'   `off_axis` and `internal_azimuth` are in radians. See
#'   [intersect_fov()] for the contour requirements.
#'
#' @return data frame containing the selected rows of `directions`, with the
#'   same columns, row names, order, and repeated observations. An empty
#'   selection retains the table structure.
#'
#' @seealso [extract_fov()], [intersect_fov()],
#'   [clip_fov()], [fibonacci_directions()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(20, 15),
#'   normalization_radius = 25,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(40, 30), geometry)
#' fov <- extract_fov(g)
#' query_directions <- fibonacci_directions(10, max_off_axis_angle = 80)
#' query_directions <- select_directions_in_fov(query_directions, fov)
select_directions_in_fov <- function(directions, fov) {
  .check_internal_directions(directions, "directions", allow_na = TRUE)
  polygon <- .spherical_field_of_view(fov, "fov")
  selected <- .directions_in_field_of_view(directions, polygon)
  directions[selected, , drop = FALSE]
}
