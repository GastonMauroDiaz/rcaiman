#' Reproject to equidistant
#'
#' @param r \linkS4class{Raster}. Only methods for images with one or three
#'   layers have been implemented.
#' @inheritParams ootb_mblt
#' @param radius Numeric integer of length one. Radius of the reprojected
#'   hemispherical image.
#'
#' @export
#'
#' @family Lens functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- apply_thr(caim$Blue, 0.5)
#' bin_equi <- reproject_to_equidistant(bin, z, a, radius = 400)
#' bin_equi <- apply_thr(bin_equi, 0.5)
#' plot(bin_equi)
#' # use write_bin(bin, "path\file_name") to have a file ready
#' # for calculating LAI with CIMES, GLA, CAN-EYE, etc.
#' }
setGeneric("reproject_to_equidistant",
           function(r, z, a, radius = 745)
             standardGeneric("reproject_to_equidistant"))

#' @rdname reproject_to_equidistant
setMethod("reproject_to_equidistant",
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
            cart <- pracma::pol2cart(as.matrix(pol))

            p <- sp::SpatialPointsDataFrame(cart[,1:2], data.frame(x = cart[,3]))
            new_r <- relative_radius_image(radius * 2)
            extent(new_r) <- extent(-pi / 2, pi / 2, -pi / 2, pi / 2)
            fun <- function(x,...) mean(x, na.rm = TRUE)
            r <- raster::rasterize(p, new_r, "x", fun)
            extent(r) <- extent(z)
            r
          }
)

#' @rdname reproject_to_equidistant
setMethod("reproject_to_equidistant",
          signature(r = "RasterStackBrick"),
          function (r, z, a, radius)
          {
            stopifnot(class(z) == "RasterLayer")
            stopifnot(class(a) == "RasterLayer")
            stopifnot(.get_max(z) <= 90)
            compareRaster(r, z)
            compareRaster(z, a)

            red <- suppressWarnings(reproject_to_equidistant(r$Red, z, a, radius))
            green <- suppressWarnings(reproject_to_equidistant(r$Green, z, a, radius))
            blue <-  suppressWarnings(reproject_to_equidistant(r$Blue, z, a, radius))
            r <- brick(red, green, blue)
            names(r) <- c("Red", "Green", "Blue")
            r
          }
)
