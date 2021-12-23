#' Read binarized images
#'
#' Wrapper functions for \code{\link[raster]{raster}}.
#'
#' @param path One-length character vector. Path to read or a binarized image.
#'
#' @export
#'
#' @seealso \code{\link{write_bin}}
#' @family Tools functions
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' m <- !is.na(z)
#' write_bin(m, "mask")
#' m_from_disk <- read_bin("mask.tif")
#' plot(m - m_from_disk)
#' }
read_bin <- function(path) {
  r <- raster(path)
  r <- is.na(r)
  if (stats::sd(r[]) == 0) r <- raster(path)
  r
}
