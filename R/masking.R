#' Image masking
#'
#' @param r \linkS4class{Raster}. The image. Values should be normalized, see
#'   \code{\link{normalize}}.
#' @param m \linkS4class{RasterLayer}. The mask, see
#'   \code{\link{mask_hs}}.
#' @param RGB Numeric vector of length three. RGB color code. Red is the default
#'   color.
#'
#' @return \linkS4class{RasterStack}
#' @export
#' @family Tools functions
#' @seealso \code{\link{mask_hs}}
#'
#' @examples
#'  r <- read_caim()
#'  z <- zenith_image(ncol(r), lens())
#'  a <- azimuth_image(z)
#'  m <- mask_hs(z, 20, 70) & mask_hs(a, 90, 180)
#'
#'  masked_caim <-  masking(normalize(r, 0, 255), m)
#'  plotRGB(masked_caim * 255)
#' \dontrun{
#'  masked_bin <- masking(apply_thr(r$Blue, 125), m)
#'  plotRGB(masked_bin * 255)
#'  }
setGeneric("masking", function(r, m, RGB = c(1,0,0))
  standardGeneric("masking"))

.masking <- function(red, green, blue, m, RGB) {
  red[!m] <- RGB[1]
  green[!m] <- RGB[2]
  blue[!m] <- RGB[3]
  stack(red, green, blue)
}

#' @rdname masking
setMethod("masking",
          signature(r = "RasterLayer"),
          function (r, m, RGB) {
            .check_if_r_was_normalized(r)
            compareRaster(r, m)
            red = green = blue <- r
            .masking(red, green, blue, m, RGB)
          }
)

#' @rdname masking
setMethod("masking",
          signature(r = "RasterStackBrick"),
          function (r, m, RGB) {
            .check_if_r_was_normalized(r)
            stopifnot(raster::nlayers(r) == 3)
            red <- raster::subset(r, 1)
            green <- raster::subset(r, 2)
            blue <- raster::subset(r, 3)
            .masking(red, green, blue, m, RGB)
          }
)
