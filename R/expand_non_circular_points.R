#' Transform points for a non-circular image expansion
#'
#' @description
#' Translate continuous image coordinates between the original frame of a
#' non-circular image and the square frame produced by
#' [expand_non_circular()].
#'
#' @param sampling_coord `data.frame` with columns `"x"` and `"y_down"`.
#'   Frame-independent columns are preserved. Columns removed by
#'   [keep_coord()] are discarded after the translation.
#' @param z numeric [terra::SpatRaster-class] produced by [zenith_image()]. It
#'   defines the square expanded image frame.
#' @param snapped_zenith_coord integer-like numeric vector of length two.
#'   Continuous image coordinates of the snapped optical center in the original
#'   image frame, in the order `c(x, y_down)`. The coordinates identify the
#'   center of the nearest 2 x 2 cell window.
#'
#' @details
#' The center of `z` has continuous image coordinates
#' `c(ncol(z) / 2, nrow(z) / 2)`. [expand_non_circular_points()] translates the
#' points so that `snapped_zenith_coord` coincides with this center.
#' [unexpand_non_circular_points()] applies the inverse translation. The
#' original image is not required for the inverse operation.
#'
#' Both functions translate every row without discarding points outside either
#' image frame. Row order, duplicates, angular directions, and nonspatial
#' columns are preserved. Derived columns tied to the previous image frame are
#' removed with [keep_coord()].
#'
#' @return
#' - `expand_non_circular_points()` returns `sampling_coord` with `"x"` and
#'   `"y_down"` expressed in the expanded image frame.
#' - `unexpand_non_circular_points()` returns `sampling_coord` with `"x"` and
#'   `"y_down"` expressed in the original image frame.
#'
#' @export
#'
#' @seealso [expand_non_circular()], [reproject_points()]
#'
#' @examples
#' z <- zenith_image(8, equidistant_coef)
#' original_coord <- data.frame(
#'   x = c(0.5, 2.5),
#'   y_down = c(0.5, 1.5),
#'   id = c("a", "b")
#' )
#' expanded_coord <- expand_non_circular_points(
#'   original_coord,
#'   z,
#'   snapped_zenith_coord = c(3, 2)
#' )
#' unexpand_non_circular_points(
#'   expanded_coord,
#'   z,
#'   snapped_zenith_coord = c(3, 2)
#' )
expand_non_circular_points <- function(sampling_coord,
                                       z,
                                       snapped_zenith_coord) {
  translation_coord <- .expand_non_circular_translation_coord(
    z,
    snapped_zenith_coord
  )
  .translate_points_coord(sampling_coord, translation_coord)
}


#' @rdname expand_non_circular_points
#' @export
unexpand_non_circular_points <- function(sampling_coord,
                                         z,
                                         snapped_zenith_coord) {
  translation_coord <- .expand_non_circular_translation_coord(
    z,
    snapped_zenith_coord
  )
  .translate_points_coord(sampling_coord, -translation_coord)
}
