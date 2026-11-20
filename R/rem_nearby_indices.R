#' Remove nearby sampling indices
#'
#' @description
#' Select a subset of cell indices so that no retained pair is closer than
#' `min_dist` in planar or spherical space.
#'
#' @details
#' When `space = "planar"`, distances are computed from cell indices and
#' directions (see [add_directions()]) are not required. When
#' `space = "spherical"`, angular distances are computed from the columns
#' selected by `angle_columns`, as in [nn_spherical()]. Use
#' `c("zenith", "azimuth")` for external directions or
#' `c("off_axis", "internal_azimuth")` for internal directions. Both are in
#' radians, as returned by [add_directions()]. A common rotation of all
#' directions preserves their angular separations.
#'
#' Selection follows input order unless `col_id` and `prioritize` specify a
#' ranking. Each retained row excludes its nearby candidates from subsequent
#' selection. Tied ranking values are processed in input order. This greedy
#' selection does not seek the largest possible subset, and changing the
#' processing order can change which indices are retained.
#'
#' In spherical space, each retained direction defines an exclusion cap
#' with angular radius `min_dist` (see [nn_spherical()]). Candidates strictly
#' inside that cap are excluded; those on its boundary satisfy the minimum
#' separation. This cap describes angular separation, independently of the
#' image projection.
#'
#' Missing directions cause an error. At polar angles of zero or `pi`, the
#' azimuth may be `NA` because the direction is still defined. A zero
#' `min_dist` retains all rows, including repeated indices.
#'
#' @param sampling_indices `data.frame` with positive integer-like columns
#'   `col` and `row`, numbered from 1. For `space = "spherical"`, it must also
#'   contain the angular columns selected by `angle_columns`, in radians.
#'   Additional columns are preserved.
#' @param min_dist numeric vector of length one. Nonnegative minimum allowed
#'   distance between retained indices. Units: pixels for `"planar"`, degrees
#'   for `"spherical"`. Angular distances must not exceed 180 degrees.
#' @param col_id numeric or character vector of length one or `NULL`. The name
#'   or position of the column in `sampling_indices` containing the optional
#'   numeric ranking value used to prioritize retention. Values must be finite
#'   and nonmissing. With `NULL`, use input order.
#' @param prioritize character vector of length one. Either `"higher"` or `"lower"`.
#'   Required when `col_id` is supplied: process higher or lower values first.
#'   Ignored when `col_id = NULL`.
#' @param space character vector of length one. Coordinate system for distances:
#'   `"planar"` (default) or `"spherical"`.
#' @inheritParams nn_spherical
#'
#' @note
#' Cell indices, angular directions, and any raster-derived values must refer to
#' the same underlying raster geometry. This relationship is not checked.
#' `angle_columns` is ignored for planar distances.
#'
#' @return The retained rows of `sampling_indices` as a `data.frame`, with all
#'   columns and their original row order, even when selection uses a ranking.
#'
#' @export
#'
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
#' sampling_indices <- expand.grid(col = seq(1, 80, 4), row = seq(1, 60, 4))
#' sampling_indices <- rem_unselected(sampling_indices, !is.na(d$zenith))
#' sampling_indices <- add_directions(sampling_indices,
#'                                   g = g, d = d)
#'
#' retained_indices <- rem_nearby_indices(sampling_indices,
#'   min_dist = 20, space = "spherical")
#' # The same angular separation can be evaluated in the internal frame.
#' retained_internal <- rem_nearby_indices(sampling_indices,
#'   min_dist = 20, space = "spherical",
#'   angle_columns = c("off_axis", "internal_azimuth"))
#'
#' rem_nearby_indices(sampling_indices, min_dist = 10, space = "planar")
#'
#' # Prefer higher DN when neighboring candidates compete for retention.
#' candidates <- data.frame(col = c(1, 2, 10), row = c(1, 1, 10),
#'                          dn = c(100, 200, 150))
#' rem_nearby_indices(candidates, min_dist = 2,
#'                    col_id = "dn", prioritize = "higher")
rem_nearby_indices <- function(sampling_indices,
                               min_dist = 3,
                               col_id = NULL,
                               prioritize = NULL,
                               space = "planar",
                               angle_columns = c("zenith", "azimuth")) {
  .check_vector(min_dist, "numeric", 1, sign = "nonnegative")
  .assert_choice(space, c("planar", "spherical"))
  .check_columns(sampling_indices, required_columns = c("row", "col"))
  .check_vector(sampling_indices$col, "integerish", sign = "positive",
                name = "sampling_indices$col")
  .check_vector(sampling_indices$row, "integerish", sign = "positive",
                name = "sampling_indices$row")

  processing_order <- seq_len(nrow(sampling_indices))
  if (!is.null(col_id)) {
    .assert_choice(prioritize, c("higher", "lower"))
    col_id <- .resolve_required_column(col_id, sampling_indices)
    rank_values <- sampling_indices[[col_id]]
    .check_vector(rank_values, "numeric", name = "ranking values")
    processing_order <- order(if (prioritize == "higher") -rank_values else
                                rank_values, processing_order)
  }

  if (space == "spherical") {
    if (min_dist > 180) {
      stop("`min_dist` must not exceed 180 degrees.", call. = FALSE)
    }
    .check_vector(angle_columns, "character", 2)
    if (anyDuplicated(angle_columns)) {
      stop("`angle_columns` must name two different columns.", call. = FALSE)
    }
    xyz <- .spherical_neighbor_vectors(sampling_indices, angle_columns)
    if (any(!stats::complete.cases(xyz))) {
      stop("Sampling directions must not be missing.", call. = FALSE)
    }
  }

  if (nrow(sampling_indices) == 0 || min_dist == 0) return(sampling_indices)
  if (space == "planar") {

    keep <- cpp_rem_nearby_planar(
      sampling_indices$col[processing_order],
      sampling_indices$row[processing_order],
      min_dist
    )

  } else {

    keep <- cpp_rem_nearby_spherical(
      xyz$x[processing_order],
      xyz$y[processing_order],
      xyz$z[processing_order],
      .degree2radian(min_dist)
    )

  }

  sampling_indices[sort(processing_order[keep]), , drop = FALSE]
}
