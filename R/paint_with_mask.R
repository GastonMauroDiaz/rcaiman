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
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' r <- normalize_min_max(g$off_axis)
#' z <- d$zenith
#' a <- d$azimuth
#' m <- !is.na(z) & !is.na(a) &
#'   z >= 20 & z <= 70 & a >= 90 & a <= 180
#' painted <- paint_with_mask(r, m)
#' terra::plotRGB(painted, scale = 1)
#' painted <- paint_with_mask(r, m, color = "blue", where = "inside")
#' terra::plotRGB(painted, scale = 1)
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



