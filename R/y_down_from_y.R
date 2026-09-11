#' Convert and manage image location representations
#'
#' @description
#' Convert between continuous image coordinates, ordinary x-y coordinates, and
#' cell indices. Enrich data frames with derived representations and remove
#' image-frame-dependent columns when they become obsolete.
#'
#' @details
#' `rcaiman` uses three representations of image locations:
#' \itemize{
#'   \item `coord`: continuous image coordinates, expressed as `x` and
#'     `y_down`. The origin is the upper-left image corner, with `x` increasing
#'     to the right and `y_down` increasing downward.
#'   \item `xy`: continuous x-y coordinates, expressed as `x` and `y`. The
#'     origin is the lower-left image corner, with `x` increasing to the right
#'     and `y` increasing upward, as in an ordinary x-y plot.
#'   \item `indices`: cell indices, expressed as `col` and `row`, numbered
#'     from 1, with columns increasing to the right and rows increasing
#'     downward.
#' }
#'
#' Throughout `rcaiman`, names ending in `_coord`, `_xy`, and `_indices` refer
#' to these representations, respectively.
#'
#' The term `point` is reserved for continuous locations (`coord` and `xy`),
#' because a cell index identifies an area rather than a point.
#'
#' The `coord` and `xy` representations share the same horizontal coordinate.
#' For an image of height \eqn{h}, their vertical coordinates are related by:
#'
#' \deqn{y_{down} = h - y}
#'
#' and:
#'
#' \deqn{y = h - y_{down}}
#'
#' Cell indices identify raster cells rather than continuous locations.
#' [coord_from_indices()] returns the continuous image coordinates of cell
#' centers. Thus, cell `col = 1`, `row = 1` has coordinates
#' `x = 0.5`, `y_down = 0.5`.
#' [xy_from_indices()] returns the same cell centers in ordinary x-y
#' coordinates, making its output suitable for functions such as
#' [graphics::points()].
#'
#' Conversely, [indices_from_coord()] assigns each continuous image location to
#' the cell containing it. This operation discards the subcell position and is
#' therefore not the inverse of [coord_from_indices()] for arbitrary
#' coordinates.
#'
#' The `add_*()` functions preserve the source representation and append a
#' requested derived representation. [add_xy()] appends only `y` because
#' `coord` and `xy` share `x`. Existing target columns are rejected instead of
#' overwritten. Derived columns remain valid only while the continuous image
#' coordinates and their image frame remain unchanged.
#'
#' [keep_coord()] retains `coord` as the authoritative representation. It
#' removes `col`, `row`, and `y`, as well as the frame-dependent predictors
#' `relative_radius` and `image_angle` (see [add_image_polar()]). It preserves
#' angular directions such as `zen` and `azi` and all other columns. The
#' function does not test whether the removed columns are consistent with
#' `coord`.
#'
#' The `xy` representation assumes the standardized raster geometry established
#' by [set_rcaiman_geometry()], with unit resolution and extent
#' `ext(0, ncol, 0, nrow)`.
#'
#' @section Functions:
#' \describe{
#'   \item{`y_down_from_y()`}{Convert `y` to `y_down`. Return a numeric vector
#'     of the same length as `y`.}
#'   \item{`y_from_y_down()`}{Convert `y_down` to `y`. Return a numeric vector
#'     of the same length as `y_down`.}
#'   \item{`coord_from_indices()`}{Convert cell indices to continuous image
#'     coordinates at cell centers. Return a `data.frame` with columns `x` and
#'     `y_down`.}
#'   \item{`indices_from_coord()`}{Assign continuous image coordinates to cells.
#'     Return a `data.frame` with columns `col` and `row`.}
#'   \item{`xy_from_indices()`}{Convert cell indices to continuous ordinary x-y
#'     coordinates at cell centers. Return a `data.frame` with columns `x` and
#'     `y`.}
#'   \item{`add_coord()`}{Append cell-center `x` and `y_down` coordinates to a
#'     `data.frame` containing `col` and `row`.}
#'   \item{`add_indices()`}{Append containing-cell `col` and `row` indices to a
#'     `data.frame` containing `x` and `y_down`.}
#'   \item{`add_xy()`}{Append `y` to a `data.frame` containing `x` and
#'     `y_down`.}
#'   \item{`keep_coord()`}{Retain `x` and `y_down` and remove columns that
#'     depend on the previous image frame.}
#' }
#'
#' @param y numeric vector. Vertical component of continuous x-y coordinates,
#'   measured upward from the lower image border.
#' @param y_down numeric vector. Vertical component of continuous image
#'   coordinates, measured downward from the upper image border.
#' @param image_height numeric vector of length one. Image height in pixels.
#' @param col,row numeric vectors of equal length. Column and row cell indices,
#'   numbered from 1.
#' @param x numeric vector. Horizontal component of continuous image
#'   coordinates, measured rightward from the left image border.
#' @param sampling_indices `data.frame` with columns `col` and `row`. Additional
#'   columns are preserved. The output columns `x` and `y_down` must not already
#'   exist.
#' @param sampling_coord `data.frame` with columns `x` and `y_down`. See
#'   *Functions* for columns added or removed by each function. All other
#'   columns are preserved.
#'
#' @return See *Functions*.
#'
#' @name y_down_from_y
#' @rdname y_down_from_y
#' @aliases y_from_y_down coord_from_indices indices_from_coord
#' @aliases xy_from_indices add_coord add_indices add_xy keep_coord
#'
#' @export
#'
#' @examples
#' image_height <- 1400
#'
#' # Convert between coord and xy representations
#' y_down <- 601.8
#' y <- y_from_y_down(y_down, image_height)
#' y_down_from_y(y, image_height)
#'
#' # Convert between coord and indices representations
#' coord <- data.frame(x = 710.3, y_down = 601.8)
#' indices <- indices_from_coord(coord$x, coord$y_down)
#' indices
#'
#' coord_from_indices(indices$col, indices$row)
#'
#' xy_from_indices(indices$col, indices$row, image_height)
#'
#' # Enrich data frames while preserving identifiers and source representations
#' sampling_indices <- data.frame(
#'   col = c(1, 5),
#'   row = c(1, 3),
#'   id = c("a", "b")
#' )
#' sampling_coord <- add_coord(sampling_indices)
#' sampling_coord <- add_indices(sampling_coord[c("x", "y_down", "id")])
#' sampling_coord <- add_xy(sampling_coord, image_height)
#' keep_coord(sampling_coord)
y_down_from_y <- function(y, image_height) {
  image_height - y
}


#' @rdname y_down_from_y
#' @export
y_from_y_down <- function(y_down, image_height) {
  image_height - y_down
}


#' @rdname y_down_from_y
#' @export
coord_from_indices <- function(col, row) {
  data.frame(
    x = col - 0.5,
    y_down = row - 0.5
  )
}


#' @rdname y_down_from_y
#' @export
indices_from_coord <- function(x, y_down) {
  data.frame(
    col = as.integer(floor(x)) + 1L,
    row = as.integer(floor(y_down)) + 1L
  )
}


#' @rdname y_down_from_y
#' @export
xy_from_indices <- function(col, row, image_height) {
  data.frame(
    x = col - 0.5,
    y = y_from_y_down(row - 0.5, image_height)
  )
}


#' @rdname y_down_from_y
#' @export
add_coord <- function(sampling_indices) {
  .check_columns(
    sampling_indices,
    required_columns = c("col", "row"),
    reserved_columns = c("x", "y_down")
  )

  coord <- coord_from_indices(
    sampling_indices$col,
    sampling_indices$row
  )
  cbind(sampling_indices, coord)
}


#' @rdname y_down_from_y
#' @export
add_indices <- function(sampling_coord) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down"),
    reserved_columns = c("col", "row")
  )

  indices <- indices_from_coord(
    sampling_coord$x,
    sampling_coord$y_down
  )
  cbind(sampling_coord, indices)
}


#' @rdname y_down_from_y
#' @export
add_xy <- function(sampling_coord, image_height) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down"),
    reserved_columns = "y"
  )

  sampling_coord$y <- y_from_y_down(
    sampling_coord$y_down,
    image_height
  )
  sampling_coord
}


#' @rdname y_down_from_y
#' @export
keep_coord <- function(sampling_coord) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down")
  )

  frame_dependent_columns <- intersect(
    c("col", "row", "y", "relative_radius", "image_angle"),
    names(sampling_coord)
  )
  sampling_coord[frame_dependent_columns] <- NULL
  sampling_coord[
    c("x", "y_down", setdiff(names(sampling_coord), c("x", "y_down")))
  ]
}
