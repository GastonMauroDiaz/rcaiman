#' Fisheye to panoramic
#'
#' Fisheye to panoramic (cylindrical projection)
#'
#' An early version of this function was used in
#' \insertCite{Diaz2021;textual}{rcaiman}.
#'
#' @inheritParams fisheye_to_equidistant
#' @inheritParams sky_grid_segmentation
#' @inheritParams extract_feature
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' pano <- fisheye_to_pano(caim, z, a)
#' plotRGB(pano %>% normalize_minmax() %>% multiply_by(255))
#' }
fisheye_to_pano<- function(r, z, a, fun = mean, angle_width = 1) {
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(r, z)
  terra::compareGeom(r, a)
  stopifnot(length(angle_width) == 1)

  .fisheye_to_pano <- function(r) {
    g <- sky_grid_segmentation(z, a, angle_width)
    blue <- extract_feature(r, g, fun, return_raster = FALSE)
    xy <- .decode_label(as.numeric(names(blue)))
    r <- matrix(NA, ncol = max(xy$sector_ID), nrow = max(xy$ring_ID))
    r <- terra::rast(r)
    terra::crs(r) <- terra::crs(g)
    terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
    cells <- terra::cellFromXY(r, as.matrix(xy) - 0.5)
    r[cells] <- blue
    terra::flip(r) #because nadir is 0 and for raster 0 is the bottom
  }

  if (terra::nlyr(r) == 1) {
    r <- suppressWarnings(.fisheye_to_pano(r))
  } else {
    layer_names <- names(r)
    r <- Map(function(r) .fisheye_to_pano(r), as.list(r))
    r <- terra::rast(r)
    names(r) <- layer_names
  }
  r
}
