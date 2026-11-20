#' Assign sky-grid labels
#'
#' @description
#' Divide the sky into zenith rings and geographic-azimuth sectors of equal
#' angular width, and label the observed directions on the image grid.
#'
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, in that order. Both layers are used.
#' @param angle_width positive numeric vector of length one. Width of both
#'   angular intervals in degrees. Must divide 90 into a whole number of
#'   rings (and consequently 360 into a whole number of sectors). Fewer than
#'   1000 rings are supported by the label encoding.
#' @param first_ring_different logical vector of length one. If `TRUE`, do not
#'   subdivide the first ring into sectors.
#'
#' @details
#' Rings extend from zenith to the horizon, and sectors cover a full turn
#' from north toward east. The acquisition orientation is already expressed
#' in `d`; it does not change the boundaries or IDs of this sky partition.
#' Partial coverage produces partial segments, without renumbering them.
#'
#' IDs encode both components as `sectorID * 1000 + ringID`. Ring and sector
#' IDs start at 1. Intervals include their upper boundary, and the first
#' interval on each axis also includes zero. See [ring_segmentation()] and
#' [sector_segmentation()].
#'
#' A cell receives zero if either required angle is missing or its zenith
#' angle exceeds 90 deg. When `first_ring_different = TRUE`, the entire first
#' ring receives ID 1000 and does not require azimuth. This includes a
#' direction exactly at zenith, whose azimuth is undefined.
#'
#' With the default subdivision, [sky_grid_centers()] supplies the angular
#' midpoint of each grid cell in radians, independently of image coverage.
#'
#' @return Single-layer [terra::SpatRaster-class] on the same grid as `d`,
#'   with integer segment labels and zero outside the segmentation. The
#'   object carries attributes `angle_width` and `first_ring_different`.
#' @seealso [sky_grid_centers()], [beckers_segmentation()]
#' @export
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
#' seg <- sky_grid_segmentation(d, 15)
#' plot(seg)
#' centers <- sky_grid_centers(15)
sky_grid_segmentation <- function(d, angle_width,
                                  first_ring_different = FALSE) {
  .check_vector(first_ring_different, "logical", 1)
  n <- .sky_interval_count(angle_width, 90)
  if (n >= 1000) {
    stop("The sky-grid label encoding requires fewer than 1000 rings.",
         call. = FALSE)
  }
  rings <- ring_segmentation(d, angle_width)
  sectors <- sector_segmentation(d, angle_width)
  g <- sectors * 1000 + rings
  g[rings == 0 | sectors == 0] <- 0
  if (first_ring_different) g[rings == 1] <- 1000
  names(g) <- paste0("Sky grid, ", angle_width, " degrees")
  attr(g, "angle_width") <- angle_width
  attr(g, "first_ring_different") <- first_ring_different
  g
}
