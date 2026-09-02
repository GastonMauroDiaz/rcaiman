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


#' Generate spherical Fibonacci directions
#'
#' @description
#' Generate angular directions on a spherical cap using a Fibonacci lattice
#' with approximately constant angular spacing. This operation is independent
#' of any raster geometry.
#'
#' @param spacing numeric vector of length one. Approximate angular separation
#'   between directions (deg).
#' @param max_zenith_angle numeric vector of length one. Maximum zenith angle
#'   delimiting the spherical cap (deg).
#'
#' @return `data.frame` with columns `"zen"` and `"azi"`, expressed in radians.
#'
#' @export
#'
#' @examples
#' directions <- fibonacci_zen_azi(10)
#' plot(directions$azi, directions$zen)
fibonacci_zen_azi <- function(spacing, max_zenith_angle = 90) {
  .check_vector(spacing, "numeric", 1, sign = "positive")
  .check_vector(max_zenith_angle, "numeric", 1, sign = "positive")
  if (max_zenith_angle > 180) {
    stop("`max_zenith_angle` cannot be greater than 180 deg.")
  }

  xyz <- .fibonacci(spacing, max_zenith_angle)
  as.data.frame(.cartesian_to_directions(xyz))
}


#' Map Fibonacci directions to cell indices
#'
#' Return cell indices and their angular directions after mapping a spherical
#' Fibonacci lattice onto the upper hemisphere of an image. See
#' [y_down_from_y()] for spatial conventions and [add_directions()] for angular
#' directions.
#'
#' Each Fibonacci direction is assigned directly to a raster cell. Output rows
#' follow the angular lattice order, and multiple directions assigned to the
#' same cell are retained (use [rem_nearby_indices()] to remove them). The
#' returned angles are the directions represented by the selected cells in `z`
#' and `a`, not the original lattice directions.
#'
#' @param m optional logical [terra::SpatRaster-class] with one layer. A binary
#'   mask with `TRUE` for selected pixels. Only indices of selected pixels are
#'   returned.
#'
#' @inheritParams add_directions
#' @inheritParams fibonacci_zen_azi
#'
#' @return `data.frame` with columns `"col"`, `"row"` and directions in
#'   radians as in [add_directions()].
#'
#' @seealso [sky_grid_centers()], [raster_grid_in_circle()]
#'
#' @export
#'
#' @examples
#' z <- zenith_image(100, equidistant_coef)
#' a <- azimuth_image(z)
#' indices <- fibonacci_indices(z, a, 30)
#' plot(z)
#' points(xy_from_indices(indices$col, indices$row, nrow(z)))
fibonacci_indices <- function(z, a, spacing, m = NULL) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
  }

  ang <- fibonacci_zen_azi(spacing)

  col_row <- indices_from_zen_azi(
    z,
    a,
    zen = ang[, "zen"],
    azi = ang[, "azi"]
  )
  col_row <- col_row[stats::complete.cases(col_row), , drop = FALSE]

  foo <- .add_directions(col_row, z, a, in_radians = TRUE)

  foo <- foo[!is.na(foo$zen),]
  if (!is.null(m)) {
    i <- extract_dn(m, foo, use_window = FALSE)[, ncol(foo) + 1]
    foo <- foo[i,]
  }
  foo
}
