#' Remove isolated black pixels
#'
#' @description
#' Replace single black pixels (`FALSE`) that are fully surrounded by white
#' pixels (`TRUE`). Uses 8-connectivity.
#'
#' @inheritParams write_bin
#'
#' @returns Logical [terra::SpatRaster-class] of one layer.
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
#' bin <- binarize_by_region(r, ring_segmentation(z, 15), "thr_isodata")
#'
#' plot(bin)
#' plot(rem_isolated_black_pixels(bin))
#' }
rem_isolated_black_pixels <- function(bin) {
  .assert_logical_mask(bin)

  # Define a 3x3 Laplacian-like kernel
  ma <- matrix(c(1, 1, 1,
                 1, -8, 1,
                 1, 1, 1), ncol = 3)

  # Identify isolated black pixels: they have 8 white neighbors
  isolated <- terra::focal(bin, w = ma) == 8

  # Set them to white
  bin[isolated] <- 1

  # Ensure return values are logical (TRUE/FALSE)
  as.logical(bin)
}
