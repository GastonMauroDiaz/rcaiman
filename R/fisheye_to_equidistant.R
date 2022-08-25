#' Fisheye to equidistant
#'
#' Fisheye to equidistant projection (also known as polar projection).
#'
#' There is no interpolation, so \code{NA} values may be generated depending on
#' both the \code{radius} argument and how much the lens projection differs from
#' the polar one. As a rule of thumb, increase \code{radius} as long as it does
#' not produce \code{NA} values on the regions to be analyzed.
#'
#' @param r \linkS4class{SpatRaster}.
#' @inheritParams ootb_mblt
#' @param radius Numeric integer of length one. Radius of the reprojected
#'   hemispherical image (i.e., the output).
#'
#' @export
#'
#' @family Lens Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- apply_thr(caim$Blue, 0.5)
#' bin_equi <- fisheye_to_equidistant(bin, z, a, radius = 400)
#' bin_equi <- apply_thr(bin_equi, 0.5)
#' plot(bin_equi)
#' # use write_bin(bin, "path\file_name") to have a file ready
#' # for calculating LAI with CIMES, GLA, CAN-EYE, etc.
#' }
fisheye_to_equidistant <- function(r, z, a, radius = 745) {
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(r, z)
  terra::compareGeom(r, a)
  stopifnot(length(radius) == 1)

  .fisheye_to_equidistant <- function(r) {
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
    r <- suppressWarnings(.fisheye_to_equidistant(r))
  } else {
    layer_names <- names(r)
    r <- Map(function(r) .fisheye_to_equidistant(r), as.list(r))
    r <- terra::rast(r)
    names(r) <- layer_names
  }
  r
}
