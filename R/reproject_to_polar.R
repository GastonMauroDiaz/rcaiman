#' Reproject to polar
#'
#' @inheritParams ootb_mblt
#' @param radius Numeric integer of length one. Radius of the reprojected
#'   hemispherical image.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- apply_thr(caim$Blue, 0.5)
#' bin_polar <- reproject_to_polar(bin, z, a, radius = 400)
#' bin_polar <- apply_thr(bin_polar, 0.5)
#' write_bin(bin_polar, "bin") #ready for CIMES, GLA, CAN-EYE, etc.
#' }
setGeneric("reproject_to_polar",
           function(r, z, a, radius = 745)
             standardGeneric("reproject_to_polar"))

#' @rdname reproject_to_polar
setMethod("reproject_to_polar",
          signature(r = "RasterLayer"),
          function (r, z, a, radius)
          {
            stopifnot(class(z) == "RasterLayer")
            stopifnot(class(a) == "RasterLayer")
            stopifnot(.get_max(z) <= 90)
            compareRaster(r, z)
            compareRaster(z, a)

            m <- !is.na(z)
            pol <- data.frame(theta = a[m] * pi / 180 + pi / 2,
                              r = z[m] * pi / 180,
                              z = r[m])
            cart <- geometry::pol2cart(as.matrix(pol))

            p <- sp::SpatialPointsDataFrame(cart[,1:2], data.frame(x = cart[,3]))
            new_r <- relative_radius_image(radius * 2)
            extent(new_r) <- extent(-pi / 2, pi / 2, -pi / 2, pi / 2)
            fun <- function(x,...) mean(x, na.rm = TRUE)
            r <- raster::rasterize(p, new_r, "x", fun)
            extent(r) <- extent(z)
            r
          }
)

#' @rdname reproject_to_polar
setMethod("reproject_to_polar",
          signature(r = "RasterStackBrick"),
          function (r, z, a, radius)
          {
            stopifnot(class(z) == "RasterLayer")
            stopifnot(class(a) == "RasterLayer")
            stopifnot(.get_max(z) <= 90)
            compareRaster(r, z)
            compareRaster(z, a)

            red <- suppressWarnings(reproject_to_polar(r$Red, z, a, radius))
            green <- suppressWarnings(reproject_to_polar(r$Green, z, a, radius))
            blue <-  suppressWarnings(reproject_to_polar(r$Blue, z, a, radius))
            r <- brick(red, green, blue)
            names(r) <- c("Red", "Green", "Blue")
            r
          }
)
