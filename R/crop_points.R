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
#'   `c(col, row)`. Cell indices are numbered from 1.
#'
#' @details
#' The upper-left boundary of the first retained cell has original image
#' coordinates `c(col - 1, row - 1)`. This value is subtracted from every input
#' point. The function translates all rows and does not discard points outside
#' the cropped image bounds.
#'
#' @return `data.frame` equal to `sampling_coord`, with `"x"` and `"y_down"`
#'   expressed in the cropped image frame.
#'
#' @export
#'
#' @examples
#' sampling_coord <- data.frame(
#'   x = c(3.5, 5.2),
#'   y_down = c(2.5, 4.8),
#'   id = c("a", "b")
#' )
#' crop_points(sampling_coord, upper_left_indices = c(4, 3))
crop_points <- function(sampling_coord, upper_left_indices) {
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
  .check_vector(
    upper_left_indices,
    "integerish",
    2,
    sign = "positive"
  )

  offset_coord <- upper_left_indices - 1
  sampling_coord$x <- sampling_coord$x - offset_coord[1]
  sampling_coord$y_down <- sampling_coord$y_down - offset_coord[2]
  sampling_coord
}
