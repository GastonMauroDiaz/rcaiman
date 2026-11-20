#' Convert between internal directions and image coordinates
#'
#' @description
#' Convert internal viewing directions to continuous image coordinates, or
#' image coordinates to internal directions, using a declared projection model.
#'
#' @section Interface stability:
#' `r lifecycle::badge("stable")`
#' The radial model interface is stable in both functions.
#'
#' `r lifecycle::badge("experimental")`
#' The `geometry$model = "anisotropic_projection"` option in
#' `internal_directions_from_coord()` is experimental.
#'
#' @details
#' Both functions use the instrumental geometry supplied in `geometry`.
#' Acquisition orientation is not involved: directions are expressed by
#' `off_axis` and clockwise `internal_azimuth`, in radians. See
#' [set_internal_geometry()] for the models and their parameters.
#'
#' `coord_from_internal_directions()` evaluates the radial polynomial and uses the
#' internal azimuth to locate each direction relative to `center_of_vision`.
#' Existing external-angle columns, if present, are not used or
#' modified. Existing `x` and `y_down` are replaced, and derived spatial columns
#' are removed with [keep_coord()].
#'
#' `internal_directions_from_coord()` supports both radial and anisotropic models using
#' the same evaluators as [set_internal_geometry()]. The radial model measures
#' radius from `center_of_vision` and inverts the projection polynomial. The
#' anisotropic model evaluates off-axis angle and internal azimuth from the
#' relative radius and image-plane angle, without requiring [add_image_polar()]
#' as a preparation step. Existing polar-predictor columns are not used.
#' It adds or replaces `off_axis` and `internal_azimuth`, preserving the input
#' coordinates and other columns. This evaluates the model at continuous
#' positions; use [add_directions()] to read directions from raster cells.
#'
#' No positions are snapped to pixels. Use [indices_from_coord()] or
#' [add_indices()] to assign the resulting positions to cells. Multiple
#' directions can fall in the same cell; input rows are not deduplicated.
#' Coordinates are not clipped to an image frame, even when `dim` is supplied.
#' Use a field-of-view contour to select directions when needed.
#'
#' Missing input coordinates produce `NA` in both output angles. In the radial
#' model, positions beyond `max_off_axis_angle` also produce `NA`. The
#' anisotropic model retains its predictions without clipping them to a
#' hemisphere or applying the radial model's angular limit. When both modeled
#' sine and cosine are zero, internal azimuth is `NA`.
#'
#' The optical-axis direction maps to
#' `center_of_vision` even when its internal azimuth is missing. Conversely,
#' that exact position has `off_axis = 0` and `internal_azimuth = NA`. The
#' opposite direction (`off_axis = pi`) has no unique radial position and returns
#' `NA` coordinates, even when the angular limit is 180 deg.
#'
#' @param directions data frame with numeric columns `off_axis` in `[0, pi]`
#'   and `internal_azimuth` in `[0, 2 * pi)`, in radians, as returned by
#'   [fibonacci_directions()] or [add_directions()]. Missing angles may be `NA`.
#' @param coord data frame with numeric columns `x` and `y_down`, in pixels,
#'   following [y_down_from_y()]. Coordinates may be fractional or negative.
#'   Missing coordinates may be `NA`.
#' @param geometry named list following [set_internal_geometry()], with
#'   `model`, `center_of_vision`, and `parameters`. For `"radial_projection"`,
#'   also supply `normalization_radius`, `normalization_angle`, and optionally
#'   `max_off_axis_angle` (90 deg when absent). `internal_directions_from_coord()` also
#'   accepts `"anisotropic_projection"`; `coord_from_internal_directions()` remains
#'   limited to the radial model. Additional fields and list classes are accepted.
#' @param dim optional positive integer-like vector of length two,
#'   `c(width, height)`. Required by `internal_directions_from_coord()` for the
#'   anisotropic model and ignored for the radial model. It defines the image
#'   semidiagonal used to normalize radius and must correspond to the frame
#'   used to fit the model and to the supplied coordinates. It does not create
#'   a raster or assign coordinates to cells.
#'
#' @return data frame with input row order, row names, and repeated observations
#'   preserved. `coord_from_internal_directions()` returns `x` and `y_down` first,
#'   followed by retained input columns. `internal_directions_from_coord()` retains
#'   the input columns and adds or replaces the two internal-angle columns.
#'
#' @seealso [add_directions()], [internal_directions_from_external()],
#'   [external_directions_from_internal()], [reproject_points()],
#'   [calc_relative_radius()], [select_directions_in_fov()]
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
#' sampling_directions <- fibonacci_directions(5)
#' sampling_directions <- select_directions_in_fov(sampling_directions, fov)
#' sampling_coord <- coord_from_internal_directions(sampling_directions, geometry)
#' sampling_indices <- indices_from_coord(sampling_coord$x, sampling_coord$y_down)
#' # Alternatively, retain directions and coordinates alongside the indices.
#' sampling_table <- add_indices(sampling_coord)
#'
#' # Evaluate the model at the centers of the selected cells.
#' cell_coord <- coord_from_indices(sampling_indices$col, sampling_indices$row)
#' cell_directions <- internal_directions_from_coord(cell_coord, geometry)
#' head(cell_directions)
#'
#' # Evaluate an anisotropic model at continuous positions, without a raster.
#' geometry <- list(
#'   model = "anisotropic_projection",
#'   center_of_vision = c(20, 15),
#'   parameters = list(
#'     off_axis_coef = c(1.1, 0.1),
#'     sin_internal_azimuth_coef = c(0, 1),
#'     cos_internal_azimuth_coef = c(-1, 0),
#'     atan_scale = 0.6
#'   )
#' )
#' coord <- data.frame(x = c(12.3, 20, 25.7), y_down = c(8.9, 15, 22.1))
#' internal_directions_from_coord(coord, geometry, dim = c(40, 30))
coord_from_internal_directions <- function(directions, geometry) {
  .check_radial_destination(geometry)
  .check_internal_directions(directions, "directions", allow_na = TRUE)
  off_axis <- directions$off_axis
  azimuth <- directions$internal_azimuth
  radius <- geometry[["normalization_radius"]] *
    .calc_relative_radius(off_axis, geometry[["parameters"]])
  center <- geometry[["center_of_vision"]]
  x <- center[1] + radius * sin(azimuth)
  y_down <- center[2] - radius * cos(azimuth)
  tolerance <- 16 * .Machine$double.eps
  axial <- !is.na(off_axis) & abs(sin(off_axis)) <= tolerance
  on_axis <- axial & off_axis < pi/2
  x[on_axis] <- center[1]
  y_down[on_axis] <- center[2]
  limit <- if ("max_off_axis_angle" %in% names(geometry))
    geometry[["max_off_axis_angle"]] else 90
  outside <- is.na(off_axis) | off_axis > limit * pi/180 + tolerance |
    (axial & off_axis > pi/2)
  x[outside] <- y_down[outside] <- NA_real_
  directions$x <- x
  directions$y_down <- y_down
  keep_coord(directions)
}

#' @rdname coord_from_internal_directions
#' @export
internal_directions_from_coord <- function(coord, geometry, dim = NULL) {
  .check_internal_geometry(geometry)
  if (geometry[["model"]] == "anisotropic_projection") {
    if (is.null(dim)) {
      stop("Supply `dim = c(width, height)` for the anisotropic model.", call. = FALSE)
    }
    .check_vector(dim, "integerish", 2, sign = "positive")
  }
  .check_columns(coord, required_columns = c("x", "y_down"))
  .check_vector(coord$x, "numeric", allow_na = TRUE, name = "coord$x")
  .check_vector(coord$y_down, "numeric", allow_na = TRUE, name = "coord$y_down")
  limit <- if ("max_off_axis_angle" %in% names(geometry))
    geometry[["max_off_axis_angle"]] else 90
  center <- geometry[["center_of_vision"]]
  dx <- coord$x - center[1]
  dy <- coord$y_down - center[2]
  angles <- switch(geometry[["model"]],
    radial_projection = .internal_radial_directions(
      dx, dy, geometry[["parameters"]], geometry[["normalization_radius"]], limit),
    anisotropic_projection = .internal_anisotropic_directions(
      dx, dy, dim, geometry[["parameters"]])
  )
  coord$off_axis <- angles[, 1]
  coord$internal_azimuth <- angles[, 2]
  coord$internal_azimuth[!is.na(coord$internal_azimuth) &
                           coord$internal_azimuth >= 2 * pi] <- 0
  coord
}
