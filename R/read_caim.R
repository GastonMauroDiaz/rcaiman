#' Read and write canopy images
#'
#' @description
#' Read a digital image as a [terra::SpatRaster-class] or write a canopy image
#' as an 8- or 16-bit unsigned-integer GeoTIFF.
#'
#' @details
#' [read_caim()] wraps [terra::rast()], so support for image formats depends on
#' the capabilities of the `terra` package. It does not crop or resample the
#' image. Use [crop_caim()] to define a processing geometry after reading.
#'
#' [write_caim()] wraps [terra::writeRaster()]. It adds the `.tif` extension to
#' `path` if missing and writes `INT1U` data when `bit_depth = 8` or `INT2U` data
#' when `bit_depth = 16`. Files follow the standard `rcaiman` conventions and
#' can be consumed by
#' [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects).
#'
#' @section Example data:
#' The package includes example files available through
#' `system.file("external", package = "rcaiman")`.
#'
#' `sensor_data.tif` represents unprocessed sensor data extracted from a RAW
#' file, analogous to the `dcraw -D` output found in the early
#' digital hemispherical photography literature, such as
#' \insertCite{Lang2010;textual}{rcaiman}. Those files can be produced outside
#' `rcaiman`, for example with [LibRaw](https://www.libraw.org/homepage)-based
#' tools such as
#' [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects).
#'
#' Preprocessing of the example images requires explicit geometry and
#' radiometry specifications. These are also included with the package as a
#' registry snapshot (see
#' [rcaiman.registry](https://github.com/GastonMauroDiaz/rcaiman.registry)).
#'
#' \describe{
#'  \item{DSCN4248.tif and DSCN4248.tsv}{Sensor data and metadata extracted from
#'  a raw file ("DSCN4606.NEF") with rcaiman-projects. Original 16-bit TIFF file
#'  was degraded to 8-bit to save space.}
#'  \item{Nikon_Coolpix5700.FCE9_raw_poly.yaml}{Registry snapshot valid for
#'  processing "DSCN4248.tif" and with parameters for a polynomial-based
#'  radiometric correction.}
#'  \item{Nikon_Coolpix5700.FCE9_raw_tsurf.yaml}{Registry snapshot valid for
#'  processing "DSCN4248.tif" and with data for supporting a trend-surface-based
#'  radiometric correction.}
#'  \item{DSCN4500.JPG}{JPEG file obtained as described in \insertCite{Diaz2018}{rcaiman}.}
#'  \item{Nikon_Coolpix5700.FCE9_jpeg.yaml}{Registry snapshot valid for
#'  processing "DSCN4500.JPG".}
#'  \item{APC_0836.jpg}{JPEG file obtained as described in \insertCite{Diaz2023}{rcaiman}.}
#'  \item{iPhone6plus.Olloclip_jpeg.yaml}{Registry snapshot valid for
#'  processing "APC_0836.jpg".}
#'  \item{sky_points.csv and sun_disk_mark.csv}{points manually digitized over the
#'  DSCN4248.tif cropped with the geometry spec "simple_method" attached to the
#'  "exif_raw" embedded metadata signature.}
#' }
#'
#'
#' @param path character vector of length one. Source image path for
#'   [read_caim()] or destination path for [write_caim()]. When reading, include
#'   the file extension; if `NULL`, an example image is returned. When writing,
#'   the `.tif` extension is enforced.
#' @param fill_na character vector of length one. Method used to fill `NA`
#'   cells. Must be `"max"`, `"min"`, `"value"`, or `"none"`. Use `"max"` to
#'   fill with the image maximum, `"min"` to fill with the image minimum,
#'   `"value"` to fill with `fill_value`, and `"none"` to preserve `NA` values.
#' @param fill_value numeric vector of length one. Value used to fill `NA`
#'   cells when `fill_na = "value"`.
#' @param caim [terra::SpatRaster-class] to write.
#' @param bit_depth numeric vector of length one. Either `8` or `16`.
#'
#' @return
#' - `read_caim()` returns a numeric [terra::SpatRaster-class]. Layers are named
#'   `"Red"`, `"Green"`, and `"Blue"` when a three-layer image is read.
#' - `write_caim()` returns nothing and is called for its side effect.
#'
#' @note
#' `rcaiman` uses terra without geographic semantics: rasters are kept with unit
#' resolution (cell size = 1) and a standardized extent `ext(0, ncol, 0, nrow)`
#' with CRS EPSG:7589 (see [set_rcaiman_geometry()]).
#'
#' File extensions are interpreted with [tools::file_ext()].
#'
#' @export
#'
#' @seealso [crop_caim()]
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' dir(system.file("external", package = "rcaiman"))
#'
#' path <- system.file("external/DSCN4248.tif", package = "rcaiman")
#' sensor_data <- read_caim(path)
#'
#' cfa_labels <- label_cfa(sensor_data)
#' cfa_pattern <- matrix(c("Green","Magenta",
#'                         "Yellow" ,"Cyan"), byrow = TRUE, ncol = 2)
#' caim <- split_cfa(sensor_data, cfa_labels, cfa_pattern)
#' caim <- map_spectral_bands(caim, list(Red = "(Magenta + Yellow)/2",
#'                                       Green = "Green",
#'                                       Blue = "Cyan"))
#'
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_jpeg.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry, "exif_jpeg", "delta_t")
#'
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path)
#' caim <- crop_caim(
#'   caim,
#'   snapped_zenith_coord = spec$snapped_zenith_coord,
#'   horizon_radius = spec$horizon_radius
#' )
#' plot(caim)
#'
#' write_caim(caim, file.path(tempdir(), "caim"), bit_depth = 8)
#' }
read_caim <- function(path = NULL,
                      fill_na = "max",
                      fill_value = NULL) {

  .check_vector(path, "character", 1, allow_null = TRUE)
  if (is.null(path)) {
    path <- system.file("external/example.tif", package = "rcaiman")
  }
  .assert_file_exists(path)
  .assert_choice(fill_na, c("max", "min", "value", "none"))
  .check_vector(fill_value, "numeric", 1, allow_null = TRUE)

  # When the image is a JPEG or a TIFF that is not a Geo-TIFF, rast produces a
  # warning about unknown extent and flip the image.
  r <- tryCatch(terra::rast(path),
                warning = function(w) {
                  # message("'terra::rast()' flips images that do not have CRS.",
                  #         "This action was countered.")
                  terra::flip(terra::rast(path))
                }) %>%
    suppressWarnings()

  r <- set_rcaiman_geometry(r)

  if (fill_na != "none") {
    value <- switch(
      fill_na,
      max = .get_max(r),
      min = .get_min(r),
      value = fill_value
    )

    r[is.na(r)] <- value
  }

  if (terra::nlyr(r) == 3) names(r) <- c("Red", "Green", "Blue")
  r
}


#' @rdname read_caim
#' @export
write_caim <- function(caim, path, bit_depth) {
  .check_vector(path, "character", 1)

  if (!any(bit_depth == 16, bit_depth == 8)) {
    stop("bit_depth should be 8 or 16.")
  }

  caim <- set_rcaiman_geometry(caim)

  file_name <- basename(path)
  file_name <- .assert_file_extension(file_name, "tif")

  if (bit_depth == 8) {
    suppressWarnings(
      terra::writeRaster(
        caim,
        file.path(dirname(path), file_name),
        filetype = "GTiff",
        datatype = "INT1U",
        overwrite = TRUE,
        NAflag = 2^8
      )
    )
  } else {
    suppressWarnings(
      terra::writeRaster(
        caim,
        file.path(dirname(path), file_name),
        filetype = "GTiff",
        datatype = "INT2U",
        overwrite = TRUE,
        NAflag = 2^16
      )
    )
  }
}



