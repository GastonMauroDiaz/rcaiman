#' Map sky-grid centers to raster coordinates
#'
#' @description
#' Return image raster coordinates and directions (see [add_directions]) for the
#' center point of each cell in a sky grid composed of circular trapezoids of
#' equal angular resolution defined by `angle_width` such as in
#' [sky_grid_segmentation()].
#'
#' @inheritParams add_directions
#' @inheritParams fibonacci_indices
#'
#' @return `data.frame` with columns `"row"`, `"col"` and directions in
#'   radians as in [add_directions()], one per grid cell.
#'
#' @seealso [fibonacci_indices()], [raster_grid_in_circle()]
#'
#' @export
#'
#' @examples
#' z <- zenith_image(100, equidistant_coef)
#' a <- azimuth_image(z)
#' pts <- sky_grid_centers(z, a, 30)
sky_grid_centers <- function(z, a, angle_width, m = NULL) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(angle_width, "numeric", 1, sign = "positive")
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
  }

  aw_rad <- .degree2radian(angle_width)
  az <- expand.grid(seq(0 + aw_rad/2, pi*2 - aw_rad/2, aw_rad),
                    seq(0 + aw_rad/2, pi/2 - aw_rad/2, aw_rad))
  sampling_points <- row_col_from_zen_azi(z, a, az[,2], az[,1])
  foo <- data.frame(sampling_points, zen = az[,2], azi = az[,1])
  if (!is.null(m)) {
    i <- extract_dn(m, foo, use_window = FALSE)[, ncol(foo) + 1]
    foo <- foo[i,]
  }
  foo
}
