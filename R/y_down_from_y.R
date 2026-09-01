#' Convert image coordinate representations
#'
#' @description
#' Convert between continuous image coordinates, ordinary x-y coordinates, and cell
#' indices.
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
#'
#' Conversely, [indices_from_coord()] assigns each continuous image location to
#' the cell containing it. This operation discards the subcell position and is
#' therefore not the inverse of [coord_from_indices()] for arbitrary
#' coordinates.
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
#'
#' @return See *Functions*.
#'
#' @name y_down_from_y
#' @rdname y_down_from_y
#' @aliases y_from_y_down coord_from_indices indices_from_coord
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
