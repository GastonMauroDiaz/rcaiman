#' Paint with mask
#'
#' @description
#' Paint image pixels inside or outside a logical mask with a solid color.
#'
#' @param r [terra::SpatRaster-class]. The image. Values should be normalized,
#'   see [normalize_min_max()]. Only images with one or three layers are
#'   supported.
#' @param color character vector of length one or numeric vector of length
#'   three. Fill color. If character, it is converted to RGB automatically. If
#'   numeric, values must be in range \eqn{[0, 1]}.
#' @param where character vector of length one Region to paint relative to `m`.
#'   Either `"outside"` (default) or `"inside"`.
#'
#' @inheritParams compute_canopy_openness
#'
#' @return numeric [terra::SpatRaster-class] with three layers and the geometry
#'   of `r`. Equal to `r`, but with pixels in the selected region painted with
#'   `color`. Single-layer inputs are replicated to allow color painting.
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
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#'
#' rspec <- get_radiometry_spec(spec, "spectral_bands")
#' caim <- cfa_to_spectral(caim, rspec$spectral_mapping)
#'
#' caim <- normalize_min_max(caim)
#'
#' m <- select_sky_region(z, 20, 70) & select_sky_region(a, 90, 180)
#' masked_caim <-  paint_with_mask(caim, m)
#' plotRGB(masked_caim * 255)
#'
#' bin <- binarize_with_thr(caim$Blue, thr_isodata(caim$Blue[]))
#' masked_bin <- paint_with_mask(bin, m)
#' plotRGB(masked_bin * 255)
#'
#' paint_with_mask(caim, m, color = c(0.2, 0.2, 0.2))  # vector
#' paint_with_mask(caim, m, color = "blue")     # name
#' paint_with_mask(caim, m, color = "#00FF00")  # hexadecimal
#' }
paint_with_mask <- function(r, m, color = "red", where = "outside") {
  .assert_spatraster(r)
  .assert_logical_mask(m)
  .assert_same_geom(r, m)
  .assert_choice(where, c("outside", "inside"))

  # Convert character color to RGB vector
  if (is.character(color) && length(color) == 1) {
    color <- grDevices::col2rgb(color) / 255
    color <- as.numeric(color)
  }

  if (!(length(color) == 3 && is.vector(color) && is.numeric(color))) {
    stop("`color` must be character vector of length one or numeric vector of length three.")
  }

  if (terra::nlyr(r) == 1) {
    red <- green <- blue <- r
  } else if (terra::nlyr(r) == 3) {
    red   <- terra::subset(r, 1)
    green <- terra::subset(r, 2)
    blue  <- terra::subset(r, 3)
  }

  if (terra::nlyr(r) == 2 | terra::nlyr(r) > 3)
    stop("`r` must have one or three layers.")

  target <- if (where == "outside") !m else m
  red[target]   <- color[1]
  green[target] <- color[2]
  blue[target]  <- color[3]

  c(red, green, blue)
}



