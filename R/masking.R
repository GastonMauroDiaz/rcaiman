#' Image masking
#'
#' @param r \linkS4class{SpatRaster}. The image. Values should be normalized,
#'   see \code{\link{normalize}}. Only methods for images with one or three
#'   layers have been implemented.
#' @param m \linkS4class{SpatRaster}. The mask, see \code{\link{mask_hs}}.
#' @param RGB Numeric vector of length three. RGB color code. Red is the default
#'   color.
#'
#' @return An object of class \linkS4class{SpatRaster} that essentially is
#'   \code{r} with the areas delimited by \code{m} --where its pixels are equal
#'   to one-- painted in a solid color. If \code{r} is a single layer image,
#'   then the layer is triplicated to allow the use of colors.
#' @export
#' @family Tools functions
#' @seealso \code{\link{mask_hs}}
#'
#'
#' @examples
#' \dontrun{
#'  r <- read_caim()
#'  z <- zenith_image(ncol(r), lens())
#'  a <- azimuth_image(z)
#'  m <- mask_hs(z, 20, 70) & mask_hs(a, 90, 180)
#'  m <- as.logical(m)
#'
#'  masked_caim <-  masking(normalize(r, 0, 255), m)
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



