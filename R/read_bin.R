#' Read binarized images
#'
#' Wrapper functions for [terra::rast()].
#'
#' @param path Character vector of length one. Path to a binarized image.
#'
#' @export
#'
#' @return An object from class [SpatRaster-class].
#'
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' m <- !is.na(z)
#' my_file <- file.path(tempdir(), "mask.tif")
#' write_bin(m, my_file)
#' m_from_disk <- read_bin(my_file)
#' plot(m - m_from_disk)
#' }
read_bin <- function(path) {
  r <- rast(path)
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  # https://spatialreference.org/ref/sr-org/7589/
  terra::crs(r) <- "epsg:7589"
  r <- is.na(r)
  if (stats::sd(r[]) == 0) r <- rast(path)
  as.logical(r)
}
