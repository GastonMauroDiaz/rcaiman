#' Read a canopy image from a file
#'
#' Reads a born-digital image (typically RGB-JPEG or RGB-TIFF) using
#' [terra::rast()] and returns a [terra::SpatRaster-class] object.
#'
#' @details
#'
#' Internally, this is a wrapper around [terra::rast()], so support for image
#' formats depends on the capabilities of the `terra` package.
#'
#'
#'
#' This function is intended for importing color hemispherical photographs, such
#' as those obtained with digital cameras equipped with fisheye lenses. For raw
#' image files (e.g., NEF, CR2), see [read_caim_raw()]. A typical workflow makes
#' use of [read_caim_raw()], [write_caim()], and [read_caim()]
#'
#' The example image was created from a raw photograph taken with a Nikon Coolpix
#' 5700 and a FC-E9 auxiliary lens, processed with the following code:
#'
#' ```
#' zenith_colrow <- c(1290, 988)/2
#' diameter <- 756
#' z <- zenith_image(diameter, lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' caim <- read_caim_raw("DSCN4606.NEF")
#' caim <- crop_caim(caim, zenith_colrow - diameter/2, diameter, diameter)
#' caim <- correct_vignetting(caim, z, c(0.0638, -0.101))
#' caim <- c(mean(caim$Y, caim$M), caim$G, caim$C)
#' caim <- fisheye_to_equidistant(caim, z, a, m, radius = 300, k = 1)
#' write_caim(caim, "example.tif", 16)
#' ```
#'
#' @param path character vector of length one.  Path to an image file, including
#'   extension. If `NULL`, an example image is returned.
#' @param zenith_colrow numeric vector of length 2. Column and row coordinates
#'   of the optical center (zenith) in pixel units.
#' @param horizon_radius integer-like numeric vector of length one. Radius of the
#'   hemispherical image in pixel units, measured from the optical center
#'   (`zenith_colrow`) to the image location corresponding to a zenith angle of
#'   90 deg (i.e. the horizon). It must be an even integer (see
#'   [rcaiman::zenith_image]).
#'
#' @return Numeric [terra::SpatRaster-class], typically with layers named
#'   `"Red"`, `"Green"`, and `"Blue"`. If the file extension does not correspond
#'   to typical JPEG or TIFF files, names will be inferred and a
#'   warning will be issued.
#'
#' @note `rcaiman` uses terra without geographic semantics: rasters are kept with
#'   unit resolution (cell size = 1) and a standardized extent
#'   `ext(0, ncol, 0, nrow)` with CRS EPSG:7589 ([set_rcaiman_geometry()]).
#'
#' @export
#'
#' @seealso [write_caim()]
#'
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' zenith_colrow <- c(1286, 986)
#' horizon_radius <- 742
#' caim <- read_caim(path, zenith_colrow, horizon_radius)
#' plot(caim$Blue)
read_caim <- function(path = NULL, zenith_colrow = NULL, horizon_radius = NULL) {

  .check_vector(path, "character", 1, allow_null = TRUE)
  if (is.null(path)) {
    path <- system.file("external/example.tif", package = "rcaiman")
  }
  .assert_file_exists(path)
  .check_vector(zenith_colrow, "numeric", length = 2, sign = "positive", allow_null = TRUE)
  .check_vector(horizon_radius, "even_integerish", length = 1, allow_null = TRUE)

  # When the image is a JPEG or a TIFF that is not a Geo-TIFF, rast produces a
  # warning about unknown extent and flip the image.
  r <- tryCatch(terra::rast(path),
                warning = function(w) {
                  message("'terra::rast()' flips images that do not have CRS.",
                          "This action was countered.")
                  terra::flip(terra::rast(path))
                }) %>%
    suppressWarnings()

  r <- set_rcaiman_geometry(r)

  if (all(!is.null(zenith_colrow), !is.null(horizon_radius))) {

    upper_left <- zenith_colrow - horizon_radius
    height <- width <- horizon_radius * 2

    r <- tryCatch(crop_caim(r, upper_left, width, height), error = function(e) NULL)
    if (is.null(r)) {
      stop(
        "The selection is outside the picture border, review `zenith_colrow`",
        "and `horizon_radius`."
      )
    }
  }

  file_ext <- filenamer::as.filename(path)$ext %>% toupper()

  if (file_ext %in% c("JPG", "JPEG", "TIF", "TIFF")) {
    if (terra::nlyr(r) == 3) names(r) <- c("Red", "Green", "Blue")
  } else {
    if (terra::nlyr(r) == 3) {
      names(r) <- c("Red", "Green", "Blue")
      warning("Layers were named presuming an RGB file")
    }
  }
  r
}



