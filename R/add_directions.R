.add_directions <- function(sampling_indices, z, a, in_radians = TRUE) {
  sampling_indices <- extract_dn(
    c(z, a),
    sampling_indices,
    use_window = FALSE
  )
  angular_columns <- ncol(sampling_indices) - c(1, 0)
  names(sampling_indices)[angular_columns] <- c("zen", "azi")
  if (in_radians) {
    sampling_indices$zen <- .degree2radian(sampling_indices$zen)
    sampling_indices$azi <- .degree2radian(sampling_indices$azi)
  }
  sampling_indices
}

#' Add viewing directions to sampling indices
#'
#' @description
#' Add internal or external viewing directions in radians to a table of cell
#' indices.
#'
#' @details
#' Directions are read from the indexed cells without interpolation or local
#' averaging. The rasters store degrees for visual inspection; this function
#' converts the extracted angles to radians for numerical analysis.
#'
#' Supply either geometry or both. Internal and external directions can also
#' be added in successive calls. Existing columns are preserved; columns that
#' would be added must not already exist.
#'
#' The function stops when a requested direction is missing, including cells
#' outside the raster. Use [rem_unselected()] to select points before adding
#' directions. At zenith and nadir, external azimuth may be `NA` because the
#' direction is defined by its zenith angle alone. This `NA` is retained. Filtering
#' with `!is.na(d$zenith)` preserves these vertical directions.
#'
#' @note
#' Use `add_directions()` when observations refer to image cells: the raster
#' supplies the direction stored at each selected cell center. If continuous
#' coordinates are first converted to indices, different positions within the
#' same cell receive the same direction.
#'
#' Use [internal_directions_from_coord()] to evaluate the geometry model directly at
#' continuous `x` and `y_down` positions, preserving subpixel differences and
#' avoiding the construction of a raster. It returns internal directions only.
#'
#' @param sampling_indices data frame with positive integer-like columns
#'   `col` and `row`, numbered from 1. Additional columns are preserved.
#' @param g optional [terra::SpatRaster-class] from
#'   [set_internal_geometry()], with layers `off_axis` and `internal_azimuth`
#'   in degrees, in that order.
#' @param d optional [terra::SpatRaster-class] from [direction_image()], with
#'   layers `zenith` and `azimuth` in degrees, in that order. At least one of
#'   `g` and `d` must be supplied. When both are supplied,
#'   they must share the same grid.
#'
#' @return data frame containing the original rows and columns, with
#'   `off_axis` and `internal_azimuth` when internal geometry is supplied,
#'   and `zenith` and `azimuth` when external geometry is supplied.
#'   All added angles are in radians. Row order and repeated indices are
#'   preserved.
#'
#' @seealso [internal_directions_from_coord()], [rem_unselected()], [extract_dn()],
#'   [indices_from_coord()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(5, 3),
#'   normalization_radius = 3,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(10, 6), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' sampling_indices <- data.frame(col = c(1, 5, 6, 7), row = c(1, 3, 3, 4))
#' m <- !is.na(g$off_axis) & !is.na(g$internal_azimuth)
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' points <- add_directions(sampling_indices, g = g)
#' points <- add_directions(points, d = d)
#' points
add_directions <- function(sampling_indices, g = NULL, d = NULL) {
  if (is.null(g) && is.null(d)) {
    stop("Supply `g`, `d`, or both.", call. = FALSE)
  }
  sources <- list(g = g, d = d)
  columns <- list(g = c("off_axis", "internal_azimuth"),
                   d = c("zenith", "azimuth"))
  for (id in names(sources)) {
    r <- sources[[id]]
    if (is.null(r)) next
    .assert_spatraster(r, name = id)
    if (!identical(names(r), columns[[id]])) {
      stop("`", id, "` must have two layers named ",
           paste(columns[[id]], collapse = " and "), ", in that order.",
           call. = FALSE)
    }
  }
  if (!is.null(g) && !is.null(d)) {
    .assert_same_geom(g, d)
  }
  result <- sampling_indices
  for (id in names(sources)) {
    r <- sources[[id]]
    if (is.null(r)) next
    result <- extract_dn(r, result, use_window = FALSE)
    first <- result[[columns[[id]][1]]]
    second <- result[[columns[[id]][2]]]
    .check_vector(first, "numeric", allow_na = TRUE, name = columns[[id]][1])
    .check_vector(second, "numeric", allow_na = TRUE, name = columns[[id]][2])
    if (any(first < 0 | first > 180, na.rm = TRUE) ||
        any(second < 0 | second >= 360, na.rm = TRUE)) {
      stop("Extracted angles must be in [0, 180] and [0, 360) degrees.",
           call. = FALSE)
    }
    missing <- is.na(first) | is.na(second)
    if (id == "d") {
      vertical <- !is.na(first) & first %in% c(0, 180)
      missing[vertical] <- FALSE
    }
    if (any(missing)) {
      stop("Some requested directions are missing in `", id,
           "`. Select sampling points with `rem_unselected()` first.",
           call. = FALSE)
    }
    result[[columns[[id]][1]]] <- .degree2radian(first)
    result[[columns[[id]][2]]] <- .degree2radian(second)
  }
  result
}
