#' Read binarized images
#'
#' Wrapper functions for \code{\link[terra]{rast}}.
#'
#' @param path One-length character vector. Path to read or a binarized image.
#'
#' @export
#'
#' @return An object from class \linkS4class{SpatRaster}.
#'
#' @seealso \code{\link{write_bin}}
#' @family Tools functions
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
  r <- is.na(r)
  if (stats::sd(r[]) == 0) r <- rast(path)
  as.logical(r)
}
