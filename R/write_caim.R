#' Write canopy image
#'
#' Wrapper function for \code{\link[terra]{writeRaster}}.
#'
#' @param caim \linkS4class{SpatRaster}.
#' @param path Character vector of length one. Path for writing the image.
#' @param bit_depth Numeric vector of length one.
#'
#' @export
#'
#' @return No return value. Called for side effects.
#'
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize(., 0, 255)
#' write_caim(caim * 2^8, file.path(tempdir(), "test_8bit"), 8)
#' write_caim(caim * 2^16, file.path(tempdir(), "test_16bit"), 16)
#' }
write_caim <- function(caim, path, bit_depth) {

  if (!any(bit_depth == 16, bit_depth == 8)) {
    stop("bit_depth should be 8 or 16.")
  }

  terra::crs(caim) <- "epsg:7589" # https://spatialreference.org/ref/sr-org/7589/
  terra::ext(caim) <- terra::ext(0, ncol(caim), 0, nrow(caim))

  file_name <- basename(path)
  file_name <- .extension(file_name, "tif")

  if (bit_depth == 8) {
    suppressWarnings(
      terra::writeRaster(caim, file.path(dirname(path), file_name),
                         filetype = "GTiff", datatype = "INT1U",
                         overwrite = TRUE)
    )
  } else {
    suppressWarnings(
      terra::writeRaster(caim, file.path(dirname(path), file_name),
                         filetype = "GTiff", datatype = "INT2U",
                         overwrite = TRUE)
    )
  }
}
