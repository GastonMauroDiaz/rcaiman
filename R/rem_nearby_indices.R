#' Remove nearby sampling indices
#'
#' @description
#' Select a subset of cell indices so that no retained pair is closer than
#' `min_dist` in planar or spherical space.
#'
#' @details
#' When `space = "planar"`, distances are computed from cell indices and
#' directions (see [add_directions()]) are not required. When
#' `space = "spherical"`, distances are computed as in [nn_spherical()]. If
#' `col_id` is provided, rows are ranked by it before indices are selected.
#'
#' @param sampling_indices `data.frame` with at least columns `"col"` and
#'   `"row"`. If `space = "spherical"`, it must also contain `"zen"` and
#'   `"azi"` (see [add_directions()]).
#' @param min_dist numeric vector of length one. Minimum allowed distance
#'   between retained indices. Units: pixels for `"planar"`, deg for
#'   `"spherical"`.
#' @param col_id numeric or character vector of length one or `NULL`. The name
#'   or position of the column in `sampling_indices` containing the optional
#'   ranking value used to prioritize retention (higher values kept first).
#' @param prioritize character vector of length one. Either `"higher"` or `"lower"`.
#'   Criterion used to prioritize retention. Ignored when `col_id = NULL`.
#' @param space character vector of length one. Coordinate system for distances:
#'   `"planar"` (default) or `"spherical"`.
#'
#' @note
#' Cell indices, angular directions, and any raster-derived values must refer to
#' the same underlying raster geometry. This relationship is not checked.
#'
#' @return The retained rows of `sampling_indices` as a `data.frame`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, orthographic_coef)
#' a <- azimuth_image(z)
#' sampling_indices <- raster_grid_in_circle(z, a, 10)
#' retained_indices <- rem_nearby_indices(
#'   sampling_indices,
#'   min_dist = 20,
#'   space = "spherical"
#' )
#'
#' plot(!is.na(z))
#' points(xy_from_indices(
#'          sampling_indices$col,
#'          sampling_indices$row,
#'          nrow(z)
#'        ), col = 2, pch = 20, cex = 0.1)
#' points(xy_from_indices(
#'          retained_indices$col,
#'          retained_indices$row,
#'          nrow(z)
#'        ), col = 1, pch = 0)
#'
#' retained_indices <- rem_nearby_indices(
#'   sampling_indices,
#'   min_dist = 100,
#'   space = "planar"
#' )
#'
#' plot(!is.na(z))
#' points(xy_from_indices(
#'          sampling_indices$col,
#'          sampling_indices$row,
#'          nrow(z)
#'        ), col = 2, pch = 20, cex = 0.1)
#' points(xy_from_indices(
#'          retained_indices$col,
#'          retained_indices$row,
#'          nrow(z)
#'        ), col = 1, pch = 0)
#'
#' r <- cie_image(
#'   z,
#'   a,
#'   c(zen_deg = 45, azi_deg = 45),
#'   as.numeric(cie_table[12, 1:5])
#' )
#' sampling_indices <- extract_dn(r, sampling_indices, use_window = FALSE)
#' head(sampling_indices)
#'
#' retained_indices <- rem_nearby_indices(
#'   sampling_indices,
#'   min_dist = 100,
#'   col_id = 5,
#'   prioritize = "higher",
#'   space = "planar"
#' )
#'
#' points(xy_from_indices(
#'          retained_indices$col,
#'          retained_indices$row,
#'          nrow(z)
#'        ), col = 1, pch = 5)
#' }
rem_nearby_indices <- function(sampling_indices,
                               min_dist = 3,
                               col_id = NULL,
                               prioritize = NULL,
                               space = "planar") {
  .check_vector(min_dist, "numeric", 1, sign = "nonnegative")
  .assert_choice(space, c("planar", "spherical"))
  required_columns <- c("row", "col")
  if (space == "spherical") {
    required_columns <- c(required_columns, "zen", "azi")
  }
  if (!is.null(col_id)) {
    # TODO (2026-09-02): `col_id` and `prioritize` are validated but do not
    # affect ordering. Retention therefore always follows input order.
    .assert_choice(prioritize, c("higher", "lower"))
    col_id <- .resolve_required_column(col_id, sampling_indices)
    required_columns <- c(required_columns, col_id)
  }
  .check_columns(sampling_indices,
                 required_columns = required_columns,
                 reserved_columns = NULL
  )

  if (space == "planar") {

    keep <- cpp_rem_nearby_planar(
      sampling_indices$col,
      sampling_indices$row,
      min_dist
    )

  } else {

    xyz <- .directions_to_cartesian(sampling_indices)

    keep <- cpp_rem_nearby_spherical(
      xyz$x,
      xyz$y,
      xyz$z,
      .degree2radian(min_dist)
    )

  }

  sampling_indices[keep, ]
}
