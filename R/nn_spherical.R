#' Find nearest neighbors in spherical coordinates
#'
#' @description
#' Identify nearest sampling points to each query point based on angular
#' distance on the unit sphere.
#'
#' @details
#' Neighbor searches are performed in a three-dimensional Cartesian
#' representation of points on the unit sphere rather than by explicitly
#' computing great-circle distances through
#' [calc_spherical_distance()].
#'
#' Directional coordinates are transformed into unit vectors:
#'
#' \deqn{
#' x = \sin(\theta)\cos(\phi)
#' }
#'
#' \deqn{
#' y = \sin(\theta)\sin(\phi)
#' }
#'
#' \deqn{
#' z = \cos(\theta)
#' }
#'
#' Neighbor identification is then carried out using Euclidean distance in
#' three-dimensional space. Search radii expressed as angular distances are
#' internally converted to chord lengths:
#'
#' \deqn{
#' d = 2\sin\left(\frac{\theta}{2}\right)
#' }
#'
#' Although chord distance and great-circle distance are different geometric
#' quantities, they are in one-to-one monotonic correspondence on the unit
#' sphere. Consequently, radius and nearest-neighbor searches yield the same
#' point ordering and neighbor selection as angular distances while avoiding
#' explicit computation of pairwise spherical distances.
#'
#' This implementation is faster than [calc_spherical_distance()].
#'
#' The angular neighborhood is defined by two parameters:
#' \itemize{
#'   \item `k`: number of nearest neighbors used by `fun`.
#'   \item `search_radius`: maximum allowed spherical distance (deg). Provide
#'   NULL  (default) to override this criterion.
#' }
#'
#' Argument `rule` controls how neighbor availability is interpreted (ignored if
#' `search_radius = NULL`):
#' \itemize{
#'   \item `"any"`: all neighbors found within `search_radius` are passed to `fun`
#'                 (even if fewer than `k`).
#'   \item `"all"`: if fewer than `k` neighbors lie within the radius, no point is returned.
#' }
#'
#' @param sampling_points `data.frame` with columns `"zen"` and `"azi"`, as in
#'   the output of [add_directions()]. Candidate neighbors of `query_points`.
#'   Additional columns are ignored.
#' @param query_points `data.frame` with columns `"zen"` and `"azi"`, as in the
#'   output of [add_directions()]. Additional columns are ignored.
#' @param k numeric vector of length one. Must be positive and
#'   interpreted as integer. Number of nearest neighbors to retrieve.
#' @param search_radius numeric vector of length one or `NULL`. Maximum angular
#'   distance (deg) to include neighbors. If `NULL`, no distance filtering is
#'   applied.
#' @param rule character vector of length one. `"any"` or `"all"`. See *Details*.
#'
#' @return List of length equal to `nrow(query_points)`. Each element contains:
#' \describe{
#'   \item{`index`}{integer vector. Indices of selected neighbors in
#'     `sampling_points`.}
#'   \item{`distance`}{numeric vector. Angular distances (radians) to each
#'     neighbor.}
#' }
#'
#' @seealso [add_directions()]
#'
#' @export
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, equidistant_coef)
#' a <- azimuth_image(z)
#' sampling_points <- fibonacci_points(z, a, 5)
#' head(sampling_points)
#' nn <- nn_spherical(sampling_points, sampling_points, k = 30, search_radius = NULL)
#'
#' for (i in sample(seq_along(nn), 10)) {
#'   plot(sampling_points[nn[[i]]$index, c("row", "col")],
#'        xlim = c(0, 1000), ylim = c(0, 1000),
#'        col = grey(normalize_min_max(nn[[i]]$distance)),
#'        asp = 1)
#'   Sys.sleep(1)
#' }
#' }
#'
nn_spherical <- function(sampling_points,
                         query_points,
                         k = 1,
                         search_radius = NULL,
                         rule = "any",
                         parallel = TRUE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1) {

  .check_vector(k, "integerish", 1, sign = "positive")
  .check_vector(
    search_radius, "numeric", 1, sign = "positive", allow_null = TRUE
  )
  if (!is.null(search_radius)) .assert_choice(rule, c("any", "all"))
  .check_columns(sampling_points, required_columns = c("zen", "azi"), )
  .check_columns(query_points, required_columns = c("zen", "azi"), )
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
  } else {
    cores <- 1
  }

  search_radius_rad <- if (!is.null(search_radius)) {
    .degree2radian(search_radius)
  } else {
    NULL
  }

  sp <- .directions_to_cartesian(sampling_points)
  qp <- .directions_to_cartesian(query_points)

  neighbors <- nn_spherical_cpp(
    sp_x = sp$x,
    sp_y = sp$y,
    sp_z = sp$z,
    qp_x = qp$x,
    qp_y = qp$y,
    qp_z = qp$z,
    k = k,
    search_radius = search_radius_rad,
    cores = cores
  )
  if (!is.null(search_radius) && rule == "all") {
    lapply(seq_along(neighbors), function(i) {
      if (length(neighbors[[i]]$index) < k) {
        neighbors[[i]]$index <<- integer(0)
        neighbors[[i]]$distance <<- integer(0)
      }
    })
  }
  neighbors
}
