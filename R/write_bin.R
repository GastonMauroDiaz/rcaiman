#' Write binarized images
#'
#' Wrapper functions for \code{\link[raster]{writeRaster}}.
#'
#' @param bin \linkS4class{RasterLayer}.
#' @inheritParams write_caim
#'
#' @export
#'
#' @return No return value. Called for side effects.
#'
#' @seealso \code{\link{read_bin}}
#' @family Tools functions
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' m <- !is.na(z)
#' my_file <- file.path(tmpDir(), "mask")
#' write_bin(m, my_file)
#' extension(my_file) <- "tif"
#' m_from_disk <- read_bin(my_file)
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
