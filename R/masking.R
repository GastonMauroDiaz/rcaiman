#' Image masking
#'
#' @param r [SpatRaster-class]. The image. Values should be normalized,
#'   see [normalize_minmax()]. Only methods for images with one or three
#'   layers have been implemented.
#' @param m [SpatRaster-class]. A mask. For hemispherical photographs,
#'   check [select_sky_vault_region()].
#' @param RGB Numeric vector of length three. RGB color code. Red is the default
#'   color.
#'
#' @return An object of class [SpatRaster-class] that is
#'   `r` but with areas where `m` is equal to zero painted in a solid
#'   color. If `r` is a single layer image, then the layer is triplicated
#'   to allow the use of color.
#' @export
#'
#' @seealso [select_sky_vault_region()]
#'
#'
#' @examples
#' \dontrun{
#'  r <- read_caim()
#'  z <- zenith_image(ncol(r), lens())
#'  a <- azimuth_image(z)
#'  m <- select_sky_vault_region(z, 20, 70) & select_sky_vault_region(a, 90, 180)
#'
#'  masked_caim <-  masking(normalize_minmax(r), m)
#'  plotRGB(masked_caim * 255)
#'
#'  masked_bin <- masking(apply_thr(r$Blue, 125), m)
#'  plotRGB(masked_bin * 255)
#'  }
#'
masking <- function (r, m, RGB = c(1,0,0)) {
  stopifnot(class(r) == "SpatRaster")
  .was_normalized(r)
  .is_single_layer_raster(m, "m")
  .is_logic_and_NA_free(m)
  terra::compareGeom(r, m)
  stopifnot(length(RGB) == 3)

  if (terra::nlyr(r) == 1) {
    red = green = blue <- r
  } else {
    stopifnot(terra::nlyr(r) == 3)
    red <- terra::subset(r, 1)
    green <- terra::subset(r, 2)
    blue <- terra::subset(r, 3)
  }
  red[!m] <- RGB[1]
  green[!m] <- RGB[2]
  blue[!m] <- RGB[3]
  c(red, green, blue)
}



