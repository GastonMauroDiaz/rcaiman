#' Perform chessboard segmentation
#'
#' Segment a raster into square regions of equal side arranged in a
#' chessboard-like pattern.
#'
#' This function divides the extent of `r` into
#' non-overlapping square segments of the given side, producing a segmentation
#' map where each segment has a unique integer label. It can be an alternative
#' to [sky_grid_segmentation()] in special cases.
#'
#' @param r [terra::SpatRaster-class]
#' @param side Numeric vector of length one. Side (in pixels) of the square
#'   segments. Must be positive and interpreted as
#'   integer.
#'
#' @return [terra::SpatRaster-class] with one layer and integer values, where
#'   each unique value corresponds to a square-segment ID.
#'
#' @export
#'
#' @examples
#' r <- set_rcaiman_geometry(rast(nrow = 20, ncol = 60))
#' values(r) <- 1
#' seg <- chessboard(r, 20)
#' plot(seg)
#' seg <- chessboard(r, 13)
#' plot(seg)
chessboard <- function(r, side) {
  .assert_spatraster(r)
  .check_vector(side, "integerish", sign = "positive")

  x <- ncol(r)/side %>% trunc()
  y <- nrow(r)/side %>% trunc()
  .r <- terra::rast(ncol = x+1, nrow = y+1)
  terra::values(.r) <- 1:terra::ncell(.r)
  .r <- terra::disagg(.r, side)
  terra::ext(.r) <- terra::ext(0,ncol(.r),0,nrow(.r))
  .r <- terra::crop(.r, r)
  terra::crs(.r) <- terra::crs(r)
  names(.r) <- "Chessborad segmentation"
  .r
}
