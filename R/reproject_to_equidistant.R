#' Reproject to equidistant
#'
#' @param r \linkS4class{SpatRaster}. Only methods for images with one or three
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
reproject_to_equidistant <- function(r, z, a, radius = 745) {
  stopifnot(class(z) == "SpatRaster")
  stopifnot(class(a) == "SpatRaster")
  stopifnot(.get_max(z) <= 90)
  terra::compareGeom(r, z)
  terra::compareGeom(z, a)

  .reproject_to_equidistant <- function(r, z, a, radius) {
    m <- !is.na(z)
    pol <- data.frame(theta = a[m] * pi / 180 + pi / 2,
                      r = z[m] * pi / 180,
                      z = r[m])
    cart <- pracma::pol2cart(as.matrix(pol))
    p <- terra::vect(cart[,1:2],
                     atts = data.frame(x = cart[,3]),
                     crs = terra::crs(r))
    new_r <- zenith_image(radius * 2, lens())
    terra::ext(new_r) <- terra::ext(-pi / 2, pi / 2, -pi / 2, pi / 2)
    fun <- function(x,...) mean(x, na.rm = TRUE)
    r <- terra::rasterize(p, new_r, "x", fun)
    terra::ext(r) <- terra::ext(z)
    r[is.na(r)] <- 0
    r
  }

  if (terra::nlyr(r) == 1) {
    r <- suppressWarnings(.reproject_to_equidistant(r, z, a, radius))
  } else {
    red <- suppressWarnings(.reproject_to_equidistant(r$Red, z, a, radius))
    green <- suppressWarnings(.reproject_to_equidistant(r$Green, z, a, radius))
    blue <-  suppressWarnings(.reproject_to_equidistant(r$Blue, z, a, radius))
    r <- c(red, green, blue)
    names(r) <- c("Red", "Green", "Blue")
    r
  }
  r
}
