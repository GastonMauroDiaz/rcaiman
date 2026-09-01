#' Remove nearby sampling points
#'
#' @description
#' Select a subset of points so that no retained pair is closer than `min_dist`
#' in planar or spherical space.
#'
#' @details
#' When `space = "planar"`, distances are computed in image coordinates and
#' directions (see [add_directions()]) are not required. When `space =
#' "spherical"`, distances are computed as in [nn_spherical()]. If `col_id` is
#' provided, points are ranked by it and retained in descending order.
#'
#' @param sampling_points `data.frame` with  least columns `"row"` and `"col"`. If
#'   `space = "spherical"`, it must also contain `"zen"` and `"azi"` (see
#'   [add_directions()]).
#' @param min_dist numeric vector of length one. Minimum allowed distance
#'   between retained points. Units: pixels for `"planar"`, deg for `"spherical"`.
#' @param col_id numeric or character vector of length one or `NULL`. The name
#'   or position of the column in `sampling_points` containing the optional
#'   ranking value used to prioritize retention (higher values kept first).
#' @param prioritize character vector of length one. Either `"higher"` or `"lower"`.
#'   criteria used to prioritize retention. Ignored when `col_id = NULL`.
#' @param space character vector of length one. Coordinate system for distances:
#'   `"planar"` (default) or `"spherical"`.
#'
#' @note
#' It is assumed that `sampling_points` were extracted from an image with the same
#' dimensions as the `z` and `a` rasters. No checks are performed.
#'
#' @return The retained points represented as a [data.frame], same as
#'   `sampling_points`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, orthographic_coef)
#' a <- azimuth_image(z)
#' sampling_points <- raster_grid_in_circle(z, a, 10)
#' foo <- rem_nearby_points(sampling_points, min_dist = 20, space = "spherical")
#'
#' plot(!is.na(z))
#' points(sampling_points$col, nrow(z) - sampling_points$row, col = 2, pch = 20, cex = 0.1)
#' points(foo$col, nrow(z) - foo$row, col = 1, pch = 0)
#'
#' foo <- rem_nearby_points(sampling_points, min_dist = 100, space = "planar")
#'
#' plot(!is.na(z))
#' points(sampling_points$col, nrow(z) - sampling_points$row, col = 2, pch = 20, cex = 0.1)
#' points(foo$col, nrow(z) - foo$row, col = 1, pch = 0)
#'
#' r <- cie_image(z, a, c(zen_deg = 45, azi_deg = 45), as.numeric(cie_table[12, 1:5]))
#' sampling_points <- extract_dn(r, sampling_points, use_window = FALSE)
#' head(sampling_points)
#'
#' foo <- rem_nearby_points(sampling_points, min_dist = 100, space = "planar",
#'                          col_id = 5, prioritize = "higher")
#'
#' points(foo$col, nrow(z) - foo$row, col = 1, pch = 5)
#' }
rem_nearby_points <- function(sampling_points,
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
    .assert_choice(prioritize, c("higher", "lower"))
    col_id <- .resolve_required_column(col_id, sampling_points)
    required_columns <- c(required_columns, col_id)
  }
  .check_columns(sampling_points,
                 required_columns = required_columns,
                 reserved_columns = NULL
  )

  if (space == "planar") {

    keep <- cpp_rem_nearby_planar(
      sampling_points$col,
      sampling_points$row,
      min_dist
    )

  } else {

    xyz <- .directions_to_cartesian(sampling_points)

    keep <- cpp_rem_nearby_spherical(
      xyz$x,
      xyz$y,
      xyz$z,
      .degree2radian(min_dist)
    )

  }

  sampling_points[keep, ]
}
