#' Write canopy image
#'
#' Wrapper function for \code{\link[raster]{writeRaster}}.
#'
#' @param caim \linkS4class{Raster}.
#' @param path Character vector of length one. Path for writing the image.
#' @param bit_depth Numeric vector of length one.
#'
#' @export
#'
#' @seealso \code{\link{write_bin}}
#' @family Tools functions
#'
#' @examples
#' \dontrun{
#' require(magrittr)
#' caim <- read_caim() %>% normalize(., 0, 255)
#' write_caim(caim * 2^8, "test_8bit", 8)
#' write_caim(caim * 2^16, "test_16bit", 16)
#' }
write_caim <- function(caim, path, bit_depth) {

  if (!any(bit_depth == 16, bit_depth == 8)) {
    stop("bit_depth should be 8 or 16.")
  }

  projection(caim) <- NA
  extent(caim) <- extent(0, ncol(caim), 0, nrow(caim))

  file_name <- basename(path)
  extension(file_name) <- "tif"

  if (bit_depth == 8) {
    suppressWarnings(
      writeRaster(caim, file.path(dirname(path), file_name),
                  format = "GTiff", datatype = "INT1U", overwrite = TRUE)
    )
  } else {
    suppressWarnings(
      writeRaster(caim, file.path(dirname(path), file_name),
                  format = "GTiff", datatype = "INT2U", overwrite = TRUE)
    )
  }
}
