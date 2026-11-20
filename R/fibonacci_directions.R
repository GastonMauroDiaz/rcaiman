.fibonacci <- function(spacing, max_angle = 90) {
  theta_max <- max_angle * pi / 180

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

#' Generate internal Fibonacci directions
#'
#' @description
#' Generate viewing directions around the optical axis using a spherical
#' Fibonacci lattice with approximately constant angular spacing.
#'
#' @details
#' `off_axis` is measured from the optical axis. `internal_azimuth` follows
#' the clock convention of [set_internal_geometry()]: zero points toward the
#' top of the sensor in landscape position, and angles increase clockwise.
#' The distribution does not depend on an image grid or external orientation.
#'
#' Directions are distributed over a spherical cap centered on the optical
#' axis, with angular radius `max_off_axis_angle` (see [nn_spherical()]).
#' At 90 deg the cap is a hemisphere; smaller radii cover a smaller angular
#' region. The generated cap does not specify the instrument's field of view.
#'
#' The same arguments produce the same directions, in the same order.
#' Approximation concerns their angular spacing, not the definition of each
#' direction. Changing `spacing` or `max_off_axis_angle` changes the distribution;
#' a smaller spherical cap is not obtained by filtering a larger one.
#' For paired measurements, generate the directions once and reuse that table.
#'
#' @param spacing positive numeric vector of length one, at most 180.
#'   Approximate angular separation between directions in degrees.
#' @param max_off_axis_angle positive numeric vector of length one, at most
#'   180. Angular limit from the optical axis in degrees. Defaults to 90.
#'
#' @return data frame with `off_axis` and `internal_azimuth` in radians,
#'   following the table convention of [add_directions()]. The number of rows
#'   depends on the angular cap and spacing; coarse spacing can yield no rows.
#'
#' @seealso [nn_spherical()], [add_directions()]
#' @export
#' @examples
#' query_directions <- fibonacci_directions(10, max_off_axis_angle = 60)
#' plot(query_directions$internal_azimuth, query_directions$off_axis,
#'      xlab = "Internal azimuth (rad)", ylab = "Off-axis angle (rad)")
#' # Find image cells around these angular queries without snapping the queries.
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(20, 15),
#'   normalization_radius = 25,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(40, 30), geometry)
#' sampling_indices <- expand.grid(col = 1:40, row = 1:30)
#' m <- !is.na(g$off_axis) & !is.na(g$internal_azimuth)
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' sampling_directions <- add_directions(sampling_indices, g)
#' neighbors <- nn_spherical(
#'   sampling_directions,
#'   query_directions,
#'   k = NULL, search_radius = 3,
#'   angle_columns = c("off_axis", "internal_azimuth"),
#'   parallel = FALSE
#' )
#' sampling_indices[neighbors[[1]]$index, ]
fibonacci_directions <- function(spacing, max_off_axis_angle = 90) {
  .check_vector(spacing, "numeric", 1, sign = "positive")
  .check_vector(max_off_axis_angle, "numeric", 1, sign = "positive")
  if (spacing > 180 || max_off_axis_angle > 180) {
    stop("`spacing` and `max_off_axis_angle` must not exceed 180 degrees.",
         call. = FALSE)
  }
  directions <- as.data.frame(.cartesian_to_directions(
    .fibonacci(spacing, max_off_axis_angle)))
  names(directions) <- c("off_axis", "internal_azimuth")
  directions
}

