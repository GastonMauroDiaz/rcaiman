#' Distance to black pixels
#'
#' @description
#' Compute the Chebyshev distance from each `TRUE` pixel to the nearest
#' `FALSE` pixel in a binary image (often depicted in black).
#'
#' @details
#' Pixels marked `FALSE` are assigned distance `0`. For pixels marked `TRUE`,
#' the returned value is the Chebyshev distance to the nearest `FALSE` pixel,
#' using 8-neighbour connectivity. Larger values indicate pixels located deeper
#' inside connected `TRUE` regions and farther from region boundaries.
#'
#' If `bin` contains no `FALSE` pixels, distances are undefined. In this case,
#' a warning is issued and a raster filled with `NA` values is returned.
#'
#' @param bin logical [terra::SpatRaster-class] of one layer.
#'
#' @return Integer [terra::SpatRaster-class] of one layer. Pixels corresponding
#'   to `FALSE` values in `bin` are assigned `0`; pixels corresponding to
#'   `TRUE` values receive their distance to the nearest `FALSE` pixel.
#'
#' @seealso [label_components()]
#'
#' @export
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
#' thr <- apply_by_direction(r, z, a, m, spacing = 10, fov = 40)
#' plot(thr$v1)
#' plot(thr$n)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' dtb <- distance_to_black(bin)
#' plot(dtb)
#' }
distance_to_black <- function(bin) {
  .assert_logical_mask(bin)

  x <- terra::as.matrix(bin, wide = TRUE)

  d <- distance_to_black_cpp(x)

  terra::setValues(bin, d)
}
