#' Label connected components
#'
#' @description
#' Assign a unique integer label to each connected component of a binary image.
#' Components are defined as groups of adjacent `TRUE` pixels using
#' 4-neighbour connectivity.
#'
#' @param bin logical [terra::SpatRaster-class] of one layer.
#'
#' @return Integer [terra::SpatRaster-class] of one layer. Pixels belonging to
#'   black (`FALSE` in `bin`) are assigned label `0`. Connected white components
#'   (`TRUE`) are assigned positive integer labels starting from `1`.
#'
#' @note
#' The implementation uses an iterative flood-fill algorithm.
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
#' seg <- label_components(bin)
#' plot(seg)
#' }
label_components <- function(bin) {
  .assert_logical_mask(bin)

  x <- terra::as.matrix(bin, wide = TRUE)

  lbl <- label_components_cpp(x)

  terra::setValues(bin, lbl)
}
