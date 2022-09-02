chessboard <- function(caim, size) {
  is(caim, "SpatRaster")
  stopifnot(length(size) == 1)
  stopifnot(.is_integerish(size))

  x <- ncol(caim)/size %>% trunc()
  y <- nrow(caim)/size %>% trunc()
  r <- terra::rast(ncol = x+1, nrow = y+1)
  terra::values(r) <- 1:terra::ncell(r)
  r <- terra::disagg(r, size)
  terra::ext(r) <- terra::ext(0,ncol(r),0,nrow(r))
  r <- terra::crop(r, caim)
  terra::crs(r) <- terra::crs(caim)
  r
}
