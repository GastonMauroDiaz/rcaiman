#' Transform points between original and cropped image frames
#'
#' @description
#' Translate continuous image coordinates to or from the frame produced by a
#' cell-aligned crop.
#'
#' @param sampling_coord `data.frame` with columns `"x"` and `"y_down"`.
#'   Frame-independent columns are preserved. Columns removed by
#'   [keep_coord()] are discarded after the translation.
#' @inheritParams crop_caim
#'
#' @details
#' The upper-left boundary of the first retained cell has original image
#' coordinates `c(col - 1, row - 1)`. [crop_points()] subtracts this value
#' from every input point, whereas
#' [uncrop_points()] adds it. The original image is not required for the inverse
#' operation.
#'
#' Supply the same `upper_left_indices` used with [crop_caim()] or
#' [uncrop_caim()]. Grid dimensions are not needed to translate points.
#'
#' Both functions translate every row without discarding points outside either
#' image frame. Row order, duplicates, angular directions, and nonspatial
#' columns are preserved. Derived columns tied to the previous image frame are
#' removed with [keep_coord()].
#'
#' @return
#' - `crop_points()` returns `sampling_coord` with `"x"` and `"y_down"`
#'   expressed in the cropped image frame.
#' - `uncrop_points()` returns `sampling_coord` with `"x"` and `"y_down"`
#'   expressed in the original image frame.
#'
#' @seealso [read_points()], [crop_caim()], [uncrop_caim()]
#'
#' @export
#'
#' @examples
#' original_coord <- data.frame(
#'   x = c(3.5, 5.2),
#'   y_down = c(2.5, 4.8),
#'   id = c("a", "b")
#' )
#' cropped_coord <- crop_points(
#'   original_coord,
#'   upper_left_indices = c(4, 3)
#' )
#' uncrop_points(cropped_coord, upper_left_indices = c(4, 3))
crop_points <- function(sampling_coord, upper_left_indices) {
  .check_vector(upper_left_indices, "integerish", 2, sign = "positive")
  offset_coord <- upper_left_indices - 1
  .translate_points_coord(sampling_coord, -offset_coord)
}


#' @rdname crop_points
#' @export
uncrop_points <- function(sampling_coord, upper_left_indices) {
  .check_vector(upper_left_indices, "integerish", 2, sign = "positive")
  offset_coord <- upper_left_indices - 1
  .translate_points_coord(sampling_coord, offset_coord)
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
  keep_coord(sampling_coord)
}
