#' Mask hemisphere
#'
#' Given a zenith or azimuth image and angle restrictions, it produces a mask.
#'
#' @param r \linkS4class{RasterLayer}. The result of a call to
#'   \code{\link{zenith_image}} or \code{\link{azimuth_image}}
#' @param from,to angle in degrees, inclusive limits.
#'
#' @export
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' m1 <- mask_hemisphere(z, 20, 70)
#' plot(m1)
#' m2 <- mask_hemisphere(a, 330,360)
#' plot(m2)
#' plot(m1 & m2)
#' plot(m1 | m2)
#'
#' # if you want 15 degress at each side of 0
#' m1 <- mask_hemisphere(a, 0, 15)
#' m2 <- mask_hemisphere(a, 345, 360)
#' plot(m1 | m2)
#'
#' # better use this
#' plot(!is.na(z))
#' # instead of this
#' plot(mask_hemisphere(z, 0, 90))
mask_hemisphere <- function(r, from, to) {
  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(from) == "numeric")
  stopifnot(class(to) == "numeric")
  stopifnot(length(from) == 1)
  stopifnot(length(to) == 1)

  m <- is.na(r)
  r[is.na(r)] <- 0
  r[r >= from & r <= to] <- NA
  r <- is.na(r)
  r[m] <- 0
  r
}


#' Image masking
#'
#' @param r \linkS4class{Raster}. The image. Values should be normalized, see
#'   \code{\link{normalize}}.
#' @param m \linkS4class{RasterLayer}. The mask, see
#'   \code{\link{mask_hemisphere}}.
#' @param RGB Numeric vector of length three. RGB color code. Red is the default
#'   color.
#'
#' @return \linkS4class{RasterStack}
#' @export
#'
#' @examples
#'  r <- read_caim()
#'  z <- zenith_image(ncol(r), lens())
#'  a <- azimuth_image(z)
#'  m <- mask_hemisphere(z, 20, 70) & mask_hemisphere(a, 90, 180)
#'
#'  masked_caim <-  masking(normalize(r, 0, 255), m)
#'  plotRGB(masked_caim * 255)
#'
#'  masked_bin <- masking(apply_thr(r$Blue, 125), m)
#'  plotRGB(masked_bin * 255)
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


#' Read binarized images
#'
#' Wrapper functions for \code{\link[raster]{raster}}.
#'
#' @param path One-length character vector. Path to read or a binarized image.
#'
#' @export
#'
#' @seealso \code{\link{write_bin}}
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' m <- !is.na(z)
#' write_bin(m, "mask")
#' m_from_disk <- read_bin("mask.tif")
#' plot(m - m_from_disk)
#' }
read_bin <- function(path) {
  r <- raster(path)
  r <- is.na(r)
  if (stats::sd(r[]) == 0) r <- raster(path)
  r
}

#' Write binarized images
#'
#' Wrapper functions for \code{\link[raster]{writeRaster}}.
#'
#' @param bin \linkS4class{RasterLayer}.
#' @inheritParams write_caim
#'
#' @export
#'
#' @seealso \code{\link{read_bin}}
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' m <- !is.na(z)
#' write_bin(m, "mask")
#' m_from_disk <- read_bin("mask.tif")
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


