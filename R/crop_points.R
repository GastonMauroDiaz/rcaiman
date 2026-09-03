#' Transform points to a cropped image frame
#'
#' @description
#' Translate continuous image coordinates from an original image frame to the
#' frame produced by a cell-aligned crop.
#'
#' @param sampling_coord `data.frame` with columns `"x"` and `"y_down"`,
#'   representing continuous image coordinates in the original image frame.
#'   Additional columns are preserved.
#' @param upper_left_indices integer-like numeric vector of length two. Cell
#'   indices of the upper-left cell included in the crop, in the order
#'   `c(col, row)`. Cell indices are numbered from 1. Alternatively, use
#'   `snapped_zenith_coord` and `horizon_radius`.
#' @param snapped_zenith_coord integer-like numeric vector of length two.
#'   Continuous image coordinates of the snapped optical center in the original
#'   image frame, in the order `c(x, y_down)`. Use with `horizon_radius` as an
#'   alternative to `upper_left_indices`.
#' @param horizon_radius positive integer-like numeric vector of length one.
#'   Radius of the declared processing geometry in cell units. Use with
#'   `snapped_zenith_coord` as an alternative to `upper_left_indices`.
#'
#' @details
#' The upper-left boundary of the first retained cell has original image
#' coordinates `c(col - 1, row - 1)`. This value is subtracted from every input
#' point. The function translates all rows and does not discard points outside
#' the cropped image bounds.
#'
#' For a crop centered on a snapped zenith, the equivalent translation is
#' `snapped_zenith_coord - horizon_radius`.
#'
#' @return `data.frame` equal to `sampling_coord`, with `"x"` and `"y_down"`
#'   expressed in the cropped image frame.
#'
#' @export
#'
#' @seealso [uncrop_points()]
#'
#' @examples
#' sampling_coord <- data.frame(
#'   x = c(3.5, 5.2),
#'   y_down = c(2.5, 4.8),
#'   id = c("a", "b")
#' )
#' crop_points(sampling_coord, upper_left_indices = c(4, 3))
#'
#' crop_points(
#'   sampling_coord,
#'   snapped_zenith_coord = c(6, 5),
#'   horizon_radius = 2
#' )
crop_points <- function(sampling_coord,
                        upper_left_indices = NULL,
                        snapped_zenith_coord = NULL,
                        horizon_radius = NULL) {
  offset_coord <- .crop_offset_coord(
    upper_left_indices,
    snapped_zenith_coord,
    horizon_radius
  )
  .translate_points_coord(sampling_coord, -offset_coord)
}


.translate_points_coord <- function(sampling_coord, translation_coord) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down")
  )
  .check_vector(
    sampling_coord$x,
    "numeric",
    name = "sampling_coord$x"
  )
  .check_vector(
    sampling_coord$y_down,
    "numeric",
    name = "sampling_coord$y_down"
  )

  sampling_coord$x <- sampling_coord$x + translation_coord[1]
  sampling_coord$y_down <-
    sampling_coord$y_down + translation_coord[2]
  sampling_coord
}
