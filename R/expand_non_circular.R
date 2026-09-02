#' Expand a non-circular image
#'
#' @description
#' Add `NA` margins to a hemispherical photograph and align its snapped zenith
#' with the center of a zenith image. This operation is intended for
#' non-circular images.
#'
#' @param caim [terra::SpatRaster-class]. Typically the output of [read_caim()].
#' @param snapped_zenith_coord integer-like numeric vector of length two.
#'   Continuous image coordinates of the optical center in the original image
#'   frame, in the order `c(x, y_down)`. The coordinates identify the center of
#'   the nearest 2 x 2 cell window.
#'
#' @inheritParams add_directions
#'
#' @return [terra::SpatRaster-class] with the same layers and pixel values as `caim`,
#'   but with `NA` margins added to center the zenith.
#'
#' @note `rcaiman` uses terra without geographic semantics: rasters are kept with
#'   unit resolution (cell size = 1) and a standardized extent
#'   `ext(0, ncol, 0, nrow)` with CRS EPSG:7589 (see [set_rcaiman_geometry()]).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Non-circular fisheye images from a smartphone with an auxiliary Lens
#' # (also applicable to non-circular fisheye images from DSLR cameras)
#' require(rcaiman.registry)
#' path <- system.file("external/iPhone6plus.Olloclip_jpeg.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_jpeg", "simple_method_rescaled")
#'
#' path <- system.file("external/APC_0836.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(spec$horizon_radius * 2,  spec$parameters)
#' a <- azimuth_image(z)
#' caim <- expand_non_circular(caim, z, spec$snapped_zenith_coord)
#' plot(caim$Blue, col = seq(0, 1, 1/255) %>% grey())
#' m <- !is.na(caim$Red) & !is.na(z)
#' plot(m, add = TRUE, alpha = 0.3, legend = FALSE)
#' }
expand_non_circular <- function(caim, z, snapped_zenith_coord) {
  .assert_spatraster(caim)
  .assert_zenith_image(z)
  .check_vector(
    snapped_zenith_coord,
    "integerish",
    2,
    sign = "positive"
  )

  zenith_xy <- c(
    snapped_zenith_coord[1],
    y_from_y_down(snapped_zenith_coord[2], nrow(caim))
  )
  delta_x <-  zenith_xy[1] - (ncol(caim) / 2)
  delta_y <-  zenith_xy[2] - (nrow(caim) / 2)
  #In which quadrant is the zenith?
  # (-)(+)|(+)(+)
  #----------------
  # (-)(-)|(+)(-)
  #
  center <- ncol(z) / 2
  xmn <- center - ((ncol(caim)/2) + delta_x)
  xmx <- center + ((ncol(caim)/2) - delta_x)
  ymn <- center - ((nrow(caim)/2) + delta_y)
  ymx <- center + ((nrow(caim)/2) - delta_y)
  e <- terra::ext(xmn, xmx, ymn, ymx)
  r <- terra::deepcopy(caim)
  terra::ext(r) <- e

  r <- terra::extend(r, z)
  terra::ext(r) <- terra::align(terra::ext(r), z)
  r <- terra::resample(r, z)
  r
}

# # ============================
# # Restricted view canopy photo
# # ============================
#
# path <- system.file("external/APC_0020.jpg", package = "rcaiman")
# caim <- read_caim(path)
# plot(caim)
# caim <- normalize_min_max(caim)
# diameter <- calc_diameter(lens(), sqrt(nrow(caim)^2 + ncol(caim)^2)/2, 90)
# z <- zenith_image(diameter, lens())
# caim <- expand_non_circular(caim, z, c(ncol(caim)/2, nrow(caim)/2))
# m <- !is.na(caim$Red)
# a <- azimuth_image(z)
# caim[!m] <- 0
#
# # Knowing from metadata that the diagonal FOV is 60.3º, and assuming
# # central perspective projection, it is possible to build a lens function with
# # radius relative to 30.15º but declaring it as 90º, for later
# # generating the zenith image rescaling it accordingly. This is a workaround
# # since central perspective cannot include the horizon.
#
# theta <- seq(0, 30.15, 0.15) * pi/180
# px <- tan(theta)/tan(30.15 * pi/180)
# theta <- normalize_min_max(theta) * pi/2
#
# fit <- lm(px ~ poly(theta, 3, raw = TRUE) - 1)
# coef(fit) %>% round(3)
# validate_projection_function(c(0.579, -0.009, 0.03))
#
# z <- zenith_image(diameter, c(0.579, -0.009, 0.03))
# ring_segmentation(z, 10) %>% plot()
# z <- normalize_min_max(z, 0, 90) * 30.15
#
# plot(caim$Blue, col = seq(0,1,1/255) %>% grey())
# m <- !is.na(caim$Red) & !is.na(z)
# plot(m, add = TRUE, alpha = 0.3, legend = FALSE)
