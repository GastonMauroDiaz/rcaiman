#' Mask image
#'
#' Given angle restrictions, produce an image for masking fisheye photos.
#'
#'
#' @inheritParams sky_grid_segmentation
#' @param zlim Numeric vector of length two. Angles in degrees. Set the zenith
#'   angle range with inclusive limits.
#' @param alim Numeric vector of length two. Angles in degrees. Set the azimuth
#'   angle range with inclusive limits.
#'
#' @return \linkS4class{RasterLayer}
#' @export
#'
#' @family masking functions
#' @seealso \code{\link{write_bin}}
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#'
#' m <- mask_image(z, a, c(20, 70), c(90, 180))
#' plot(m)
#'
#' m1 <- mask_image(z, a, alim = c(90, 180))
#' plot(m1)
#'
#' m2 <- mask_image(z, zlim = c(20, 70))
#' plot(m2)
#'
#' plot(m1 & m2)
#'
#' m <- mask_image(z)
#' plot(m)
mask_image <- function(z,
                       a = NULL,
                       zlim = NULL,
                       alim = NULL) {
  no_data_area <- is.na(z)

  if (all(is.null(zlim), is.null(alim))) {
    m <- !is.na(z)
  } else {
    if (all(!is.null(zlim), !is.null(alim))) {
      stopifnot(length(zlim) == 2)
      stopifnot(length(alim) == 2)

      stopifnot(all(zlim[1] >= 0, zlim[2] <= 90))
      stopifnot(all(alim[1] >= 0, alim[2] <= 360))

      z[is.na(z)] <- 0
      a[is.na(a)] <- 0
      z[z >= zlim[1] & z <= zlim[2]] <- NA
      a[a >= alim[1] & a <= alim[2]] <- NA
      m <- is.na(z) + is.na(a)
      m[m == 2] <- NA
      m <- is.na(m)
    } else {
      if (!is.null(zlim)) {
        stopifnot(length(zlim) == 2)
        stopifnot(all(zlim[1] >= 0, zlim[2] <= 90))

        z[is.na(z)] <- 0

        z[z >= zlim[1] & z <= zlim[2]] <- NA
        m <- is.na(z)
      } else {
        stopifnot(length(alim) == 2)
        stopifnot(all(alim[1] >= 0, alim[2] <= 360))

        a[is.na(a)] <- 0
        a[a >= alim[1] & a <= alim[2]] <- NA
        m <- is.na(a)
      }
    }
  }
  # fix inclusion of the area outside the circle if zmin is 0
  m[no_data_area] <- 0
  m
}


#' Image masking
#'
#' @param r \linkS4class{Raster}. The image. Values should be normalized, see
#'   \code{\link{normalize}}.
#' @param m \linkS4class{RasterLayer}. The mask, see \code{\link{mask_image}}.
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
#'  m <- mask_image(z, a, c(20, 70), c(90, 180))
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
#' m <- mask_image(z)
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
#' m <- mask_image(z)
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


