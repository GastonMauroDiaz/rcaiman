#' Write canopy image
#'
#' @description
#' Wrapper around [terra::writeRaster()] that writes a canopy image as GeoTIFF
#' with 8- or 16-bit unsigned integers, setting CRS and extent.
#'
#' @details
#' Adds the `.tif` extension to `path` if missing. The CRS is set to EPSG:7589
#' and the extent to `[0, ncol] × [0, nrow]` in pixel units
#' (see [set_rcaiman_geometry()]). Data are written as `INT1U` when `bit_depth = 8`
#' and `INT2U` when `bit_depth = 16`.
#'
#' Files produced by this function follow the standard rcaiman file conventions.
#' These conventions enable seamless integration with
#' [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects).
#'
#' @note
#' File extensions are interpreted using [tools::file_ext()]. If `path` has no
#' extension, the default extension `"tif"` is appended automatically.
#'
#' @param caim [terra::SpatRaster-class].
#' @param path character vector of length one. Destination file path (extension
#'   `.tif` will be enforced).
#' @param bit_depth numeric vector of length one. Either `8` or `16`.
#'
#' @return No return value. Called for side effects.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' # The last bin must be keep unuse otherwise it is turned black
#' caim <- normalize_min_max(caim) * (2^8 - 1)
#' caim <- round(caim)
#' caim
#' write_caim(caim * 2^8,  file.path(tempdir(), "test_8bit"),  8)
#' }
write_caim <- function(caim, path, bit_depth) {
  .check_vector(path, "character", 1)

  if (!any(bit_depth == 16, bit_depth == 8)) {
    stop("bit_depth should be 8 or 16.")
  }

  caim <- set_rcaiman_geometry(caim)

  file_name <- basename(path)
  file_name <-  .assert_file_extension(file_name, "tif")

  if (bit_depth == 8) {
    suppressWarnings(
      terra::writeRaster(caim, file.path(dirname(path), file_name),
                         filetype = "GTiff", datatype = "INT1U",
                         overwrite = TRUE, NAflag = 2^8)
    )
  } else {
    suppressWarnings(
      terra::writeRaster(caim, file.path(dirname(path), file_name),
                         filetype = "GTiff", datatype = "INT2U",
                         overwrite = TRUE, NAflag = 2^16)
    )
  }
}
