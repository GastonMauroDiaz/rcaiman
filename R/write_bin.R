#' Write binarized images
#'
#' Wrapper functions for \code{\link[raster]{writeRaster}}.
#'
#' @param bin \linkS4class{RasterLayer}.
#' @inheritParams write_caim
#'
#' @export
#'
#' @seealso \code{\link{read_bin}}
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
write_bin <- function(bin, path) {
  stopifnot(max(bin[], na.rm = TRUE) <= 1)

  file_name <- basename(path)
  extension(file_name) <- "tif"

  projection(bin) <- NA
  extent(bin) <- extent(0, ncol(bin), 0, nrow(bin))

  suppressWarnings(
    writeRaster(bin * 255, file.path(dirname(path), file_name),
                format = "GTiff", datatype = "INT1U", overwrite = TRUE)
  )
}
