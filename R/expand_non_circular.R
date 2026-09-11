#' Place a non-circular image in a square frame
#'
#' @description
#' Translate a hemispherical photograph so that its snapped zenith coincides
#' with the center of a zenith image, and return exactly the square processing
#' frame defined by that zenith image.
#'
#' @details
#' The photograph is translated along cell boundaries. Parts of the square
#' processing frame without source-image coverage are filled with `NA`, while
#' parts of the source image outside that frame are discarded. Thus, depending
#' on the input dimensions and zenith position, the operation can require both
#' padding and cropping.
#'
#' Since `z` has even dimensions and `snapped_zenith_coord` identifies an
#' integer grid intersection, the translation is an integer number of cells.
#' The implementation extends and crops the translated raster without
#' resampling or interpolating its values.
#'
#' @param caim [terra::SpatRaster-class]. Typically the output of [read_caim()].
#' @param snapped_zenith_coord integer-like numeric vector of length two.
#'   Continuous image coordinates of the optical center in the original image
#'   frame, in the order `c(x, y_down)`. The coordinates identify the center of
#'   the nearest \eqn{2 \times 2}  cell window.
#'
#' @inheritParams add_directions
#'
#' @return [terra::SpatRaster-class] with the same geometry as `z` and the same
#'   layers as `caim`. Covered cells retain their original values; uncovered
#'   cells are `NA`.
#'
#' @note `rcaiman` uses terra without geographic semantics: rasters are kept with
#'   unit resolution (cell size = 1) and a standardized extent
#'   `ext(0, ncol, 0, nrow)` with CRS EPSG:7589 (see [set_rcaiman_geometry()]).
#'
#' @seealso [expand_non_circular_points()]
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
  translation_coord <- .expand_non_circular_translation_coord(
    z,
    snapped_zenith_coord
  )

  xmn <- translation_coord[1]
  xmx <- xmn + ncol(caim)
  ymx <- nrow(z) - translation_coord[2]
  ymn <- ymx - nrow(caim)
  e <- terra::ext(xmn, xmx, ymn, ymx)
  r <- terra::deepcopy(caim)
  terra::ext(r) <- e

  r <- terra::extend(r, z)
  # Correction (2026-09-04): extend() returns the union of both extents. Crop
  # that union to make z the exact processing frame. The grids are aligned by
  # the integer translation above, so no resampling is required.
  r <- terra::crop(r, z)
  r
}


.expand_non_circular_translation_coord <- function(z,
                                                    snapped_zenith_coord) {
  .assert_zenith_image(z)
  .check_vector(
    snapped_zenith_coord,
    "integerish",
    2,
    sign = "positive"
  )

  c(ncol(z) / 2, nrow(z) / 2) - snapped_zenith_coord
}
