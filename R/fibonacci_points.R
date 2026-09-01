.fibonacci <- function(spacing, max_zenith_angle = 90) {
  theta_max <- max_zenith_angle * pi / 180

  # Convert target cap area
  cap_area <- 2 * pi * (1 - cos(theta_max))

  theta <- spacing * pi / 180
  d <- 2 * sin(theta / 2)

  # Estimate number of points based on area
  N <- round(cap_area / (d^2))

  phi <- (1 + sqrt(5)) / 2

  pts <- matrix(NA_real_, nrow = N, ncol = 3)

  for (k in seq_len(N)) {
    # Uniform sampling over spherical cap
    zc <- 1 - (k - 0.5) / N * (1 - cos(theta_max))

    r <- sqrt(1 - zc * zc)
    az <- 2 * pi * ((k - 1) / phi %% 1)

    x <- r * cos(az)
    y <- r * sin(az)

    pts[k, ] <- c(x, y, zc)
  }

  pts
}


#' Map Fibonacci points to continuous image coordinates
#'
#' Return continuous image coordinates (see [xy_from_col_row()]) and directions (see
#' [add_directions]) for points distributed on the upper hemisphere using a
#' spherical Fibonacci lattice with approximately constant angular spacing.
#'
#' @param spacing numeric vector of length one. Angular separation (deg) between
#'   points.
#' @param m optional logical [terra::SpatRaster-class] with one layer. A binary
#'   mask with `TRUE` for selected pixels. Only points over a selected pixel are
#'   returned.
#'
#' @inheritParams add_directions
#'
#' @return `data.frame` with columns `"row"`, `"col"` and directions in
#'   radians as in [add_directions()].
#'
#' @seealso [sky_grid_centers()], [raster_grid_in_circle()]
#'
#' @export
#'
#' @examples
#' z <- zenith_image(100, equidistant_coef)
#' a <- azimuth_image(z)
#' pts <- fibonacci_points(z, a, 30)
#' \dontrun{
#' display_caim(is.na(z), sampling_points = pts)
#' }
fibonacci_points <- function(z, a, spacing, m = NULL) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(spacing, "numeric", 1, sign = "positive")
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
  }

  # generate Fibonacci points on hemisphere
  pts <- .fibonacci(spacing)
  # translate to directions
  ang <- .cartesian_to_directions(pts)

  row_col <- row_col_from_zen_azi(z, a, zen = ang[, "zen"], azi = ang[, "azi"])
  # data.frame(row_col, zen = ang[i, "zen"], azi = ang[i, "azi"]) # this does not work
  foo <- .add_directions(row_col, z, a, in_radians = TRUE) # this works
  foo <- foo[!is.na(foo$zen),]
  if (!is.null(m)) {
    i <- extract_dn(m, foo, use_window = FALSE)[, ncol(foo) + 1]
    foo <- foo[i,]
  }
  foo
}
