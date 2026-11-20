#' Label connected components
#'
#' @description
#' Assign a unique integer label to each connected component of a binary image.
#' Components are defined as groups of adjacent `TRUE` pixels using
#' 4-neighbour connectivity.
#'
#' @param bin logical [terra::SpatRaster-class] of one layer.
#'
#' @return Integer [terra::SpatRaster-class] of one layer. Pixels belonging to
#'   black (`FALSE` in `bin`) are assigned label `0`. Connected white components
#'   (`TRUE`) are assigned positive integer labels starting from `1`.
#'
#' @note
#' The implementation uses an iterative flood-fill algorithm.
#'
#' @export
#'
#' @examples
#' # Two square patches touch only at a corner; a third patch is separate.
#' bin <- matrix(c(
#'   0, 0, 0, 0, 0, 0, 0, 0,
#'   0, 1, 1, 0, 0, 0, 1, 0,
#'   0, 1, 1, 0, 0, 0, 1, 0,
#'   0, 0, 0, 1, 1, 0, 1, 0,
#'   0, 0, 0, 1, 1, 0, 0, 0,
#'   0, 0, 0, 0, 0, 0, 0, 0
#' ), nrow = 6, byrow = TRUE)
#' bin <- set_rcaiman_geometry(rast(bin)) > 0
#'
#' seg <- label_components(bin)
#' plot(c(bin, seg), main = c("Binary image", "Connected components"))
#'
#' # Four-neighbour connectivity keeps diagonally touching patches separate.
#' # Labels 1--3 identify the components; label 0 is the background.
#' table(seg[])
label_components <- function(bin) {
  .assert_logical_mask(bin)

  x <- terra::as.matrix(bin, wide = TRUE)

  lbl <- label_components_cpp(x)

  terra::setValues(bin, lbl)
}
