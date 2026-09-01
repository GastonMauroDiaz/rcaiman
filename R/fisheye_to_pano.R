#' Fisheye to panoramic
#'
#' @description
#' Reprojects a fisheye view into a panoramic view using a
#' cylindrical projection. The output is standardized so that rows correspond
#' to zenith angle bands and columns to azimuthal sectors.
#'
#' @details
#' This function computes a cylindrical projection by aggregating pixel values
#' according to their zenith and azimuth angles. Internally, it creates a
#' segmentation grid with [sky_grid_segmentation()] and applies [compute_mean()].
#'
#' @note
#' An early version of this function was used in
#' \insertCite{Diaz2021;textual}{rcaiman}.
#'
#' @inheritParams reproject_image
#' @inheritParams add_directions
#'
#' @returns
#' [terra::SpatRaster-class] with rows representing zenith angle bands and
#' columns representing azimuthal sectors. The number of layers and names
#' matches that of the input `r`.
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @seealso [panoramic()]
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' pano <- fisheye_to_pano(caim, z, a)
#' plot(pano$Cyan %>% normalize_min_max() %>% multiply_by(255))
#' }
fisheye_to_pano <- function(r, z, a, angle_width = 1) {
  .assert_spatraster(r, squared = TRUE)
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .assert_same_geom(z, r)

  .check_vector(angle_width, "numeric", 1, sign = "positive")

  .fisheye_to_pano <- function(r) {
    g <- sky_grid_segmentation(z, a, angle_width)
    blue <- computer_mean(r, g, return = "vector")
    xy <- .decode_label(as.numeric(names(blue)))
    r <- matrix(NA, ncol = max(xy$sector_ID), nrow = max(xy$ring_ID))
    r <- terra::rast(r)
    terra::crs(r) <- terra::crs(g)
    terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
    cells <- terra::cellFromXY(r, as.matrix(xy) - 0.5)
    r[cells] <- blue
    terra::flip(r) #because nadir is 0 and in rasters 0 is the bottom
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
