#' Transform points to an original image frame
#'
#' @description
#' Translate continuous image coordinates from a cropped image frame to the
#' corresponding original image frame.
#'
#' @param sampling_coord `data.frame` with columns `"x"` and `"y_down"`,
#'   representing continuous image coordinates in the cropped image frame.
#'   Additional columns are preserved.
#' @inheritParams crop_points
#'
#' @details
#' This function is the inverse of [crop_points()]. It adds the continuous image
#' coordinates of the crop origin to every input point. The original image is
#' not required because the translation is fully determined by either
#' `upper_left_indices` or `snapped_zenith_coord` and `horizon_radius`.
#'
#' The function translates all rows and does not discard points outside the
#' original image bounds.
#'
#' @return `data.frame` equal to `sampling_coord`, with `"x"` and `"y_down"`
#'   expressed in the original image frame.
#'
#' @export
#'
#' @examples
#' cropped_coord <- data.frame(
#'   x = c(0.5, 2.2),
#'   y_down = c(0.5, 2.8),
#'   id = c("a", "b")
#' )
#' uncrop_points(cropped_coord, upper_left_indices = c(4, 3))
uncrop_points <- function(sampling_coord,
                          upper_left_indices = NULL,
                          snapped_zenith_coord = NULL,
                          horizon_radius = NULL) {
  offset_coord <- .crop_offset_coord(
    upper_left_indices,
    snapped_zenith_coord,
    horizon_radius
  )
  .translate_points_coord(sampling_coord, offset_coord)
}
