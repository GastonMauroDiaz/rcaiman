#' Read binarized images
#'
#' Wrapper functions for \code{\link[raster]{raster}}.
#'
#' @param path One-length character vector. Path to read or a binarized image.
#'
#' @export
#'
#' @return An object from class \linkS4class{RasterLayer}.
#'
#' @seealso \code{\link{write_bin}}
#' @family Tools functions
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' m <- !is.na(z)
#' my_file <- file.path(tmpDir(), "mask.tif")
#' write_bin(m, my_file)
#' m_from_disk <- read_bin(my_file)
#' plot(m - m_from_disk)
#' }
read_bin <- function(path) {
  r <- raster(path)
  r <- is.na(r)
  if (stats::sd(r[]) == 0) r <- raster(path)
  r
}
