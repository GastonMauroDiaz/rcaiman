#' Write binarized images
#'
#' Wrapper functions for \code{\link[terra]{writeRaster}}.
#'
#' @param bin \linkS4class{SpatRaster}.
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
#' my_file <- file.path(tempdir(), "mask")
#' write_bin(m, my_file)
#' my_file <- as.filename(my_file) %>%
#'   insert(., ext = "tif", replace = TRUE) %>%
#'   as.character()
#' m_from_disk <- read_bin(my_file)
#' plot(m - m_from_disk)
#' }
write_bin <- function(bin, path) {
  stopifnot(.get_max(bin) <= 1)

  file_name <- basename(path)
  file_name <-  .extension(file_name, "tif")

  terra::crs(bin) <- "epsg:7589" # https://spatialreference.org/ref/sr-org/7589/
  terra::ext(bin) <- terra::ext(0, ncol(bin), 0, nrow(bin))

  suppressWarnings(
    writeRaster(bin * 255, file.path(dirname(path), file_name),
                filetype = "GTiff", datatype = "INT1U", overwrite = TRUE)
  )
}
