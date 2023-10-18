#' Chessboard segmentation
#'
#' @inheritParams polar_qtree
#' @param size Numeric vector of length one. Size of the square segments.
#'
#' @return A single layer image of the class [SpatRaster-class] with integer
#'   values.
#' @export
#'
#' @family Segmentation Functions
#'
#' @examples
#' caim <- read_caim()
#' seg <- chessboard(caim, 20)
#' plot(caim$Blue)
#' plot(extract_feature(caim$Blue, seg))
chessboard <- function(r, size) {
  is(r, "SpatRaster")
  stopifnot(length(size) == 1)
  stopifnot(.is_integerish(size))

  x <- ncol(r)/size %>% trunc()
  y <- nrow(r)/size %>% trunc()
  .r <- terra::rast(ncol = x+1, nrow = y+1)
  terra::values(.r) <- 1:terra::ncell(.r)
  .r <- terra::disagg(.r, size)
  terra::ext(.r) <- terra::ext(0,ncol(.r),0,nrow(.r))
  .r <- terra::crop(.r, r)
  terra::crs(.r) <- terra::crs(r)
  .r
}
