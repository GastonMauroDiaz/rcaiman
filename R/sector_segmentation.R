#' Assign azimuth-sector labels
#'
#' @description
#' Divide the sky into sectors of equal geographic-azimuth width and label
#' the observed directions on the image grid.
#'
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, in that order. Only `azimuth` is used;
#'   `zenith` is ignored, including its missing values.
#' @param angle_width positive numeric vector of length one. Sector width in
#'   degrees. Must divide 360 into a whole number of sectors.
#'
#' @details
#' Azimuth is measured from north toward east, as in [direction_image()].
#' Sector boundaries and IDs are independent of the acquisition orientation
#' and the fraction of sky covered by the image.
#'
#' Intervals include their upper boundary; the first also includes zero.
#' With a width of 30 deg, azimuths from 0 through 30 deg belong to sector 1,
#' and angles greater than 30 through 60 deg belong to sector 2.
#' Azimuths must be in `[0, 360)` degrees.
#'
#' Missing azimuths receive zero, including the undefined azimuth at zenith
#' or nadir. Because zenith angles are ignored, this function also assigns
#' sectors below the horizon when their azimuth is available. Use
#' [sky_grid_segmentation()] to restrict the partition to the sky hemisphere.
#'
#' @return Single-layer [terra::SpatRaster-class] on the same grid as `d`,
#'   with sector IDs starting at 1. Zero denotes cells without an azimuth.
#' @seealso [ring_segmentation()], [sky_grid_segmentation()]
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
#' sectors <- sector_segmentation(d, 30)
#' plot(sectors)
sector_segmentation <- function(d, angle_width) {
  a <- .sky_segmentation_layer(d, "azimuth")
  n <- .sky_interval_count(angle_width, 360)
  bounds <- seq(0, 360, length.out = n + 1)
  sectors <- terra::classify(a,
    cbind(head(bounds, -1), tail(bounds, -1), seq_len(n)),
    right = TRUE, include.lowest = TRUE, others = 0)
  sectors[is.na(sectors)] <- 0
  names(sectors) <- "Sectors"
  sectors
}
