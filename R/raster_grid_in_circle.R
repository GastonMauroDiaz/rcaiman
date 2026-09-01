#' Generate raster grid points within a circular mask
#'
#' @description
#' Return continuous inage coordinates (see [xy_from_col_row()]) for points
#' regularly spaced on the image plane domain and restricted to the interior of
#' the inscribed circle.
#'
#' @details
#' This function operates purely in image space. Points are generated every
#' `spacing` pixels along rows and columns, and only those strictly inside the
#' inscribed circle are retained.
#'
#' Unlike [fibonacci_points()] and [sky_grid_centers()], which distribute
#' points in angular space and then project them to the raster, this function
#' performs sampling directly on the raster grid.
#'
#' @param spacing numeric vector of length one. Must be positive and interpreted
#'   as integer. Distance in pixels between sampled points.
#'
#' @inheritParams add_directions
#' @inheritParams fibonacci_points
#'
#' @return `data.frame` with columns `"row"`, `"col"` and directions in
#'   radians as in [add_directions()].
#'
#' @seealso [fibonacci_points()], [sky_grid_centers()]
#'
#' @export
#'
#' @examples
#' z <- zenith_image(100, equidistant_coef)
#' a <- azimuth_image(z)
#' pts <- raster_grid_in_circle(z, a, 10)
#' head(pts)
raster_grid_in_circle <- function(z, a, spacing, m = NULL) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
  }

  side <- ncol(z)
  coords <- seq(1, side, by = spacing)

  foo <- expand.grid(row = coords, col = coords)

  c0 <- r0 <- side / 2

  r <- sqrt((foo$row - c0)^2 + (foo$col - c0)^2)

  foo <- .add_directions(foo[r < r0, ], z, a, in_radians = TRUE)
  foo <- foo[!is.na(foo$zen),]
  if (!is.null(m)) {
    i <- extract_dn(m, foo, use_window = FALSE)[, ncol(foo) + 1]
    foo <- foo[i,]
  }
  foo
}
