#' Find nearest neighbors in spherical coordinates
#'
#' @description
#' Select sampling directions around each query direction by angular distance.
#'
#' @details
#' Directions supplied by [add_directions()] are in radians. Select their
#' columns with `angle_columns`: `c("zenith", "azimuth")` for external
#' directions or `c("off_axis", "internal_azimuth")` for internal directions.
#' Both tables must use the same reference frame. A common rotation of all
#' directions preserves angular distances and neighbor selection.
#'
#' A spherical cap contains all directions within a given angular distance
#' of a central direction. Here, each query defines the center and a supplied
#' `search_radius` defines the angular radius of the cap. A radius of 90 deg
#' defines a hemisphere. Image coverage and masks used to select the input
#' samples determine which observations are available within the cap.
#'
#' With a numeric `k`, return up to `k` nearest neighbors, restricted by
#' `search_radius` when supplied. Set `k = NULL` and supply a radius to
#' retrieve all points within that angular neighborhood. This performs a
#' radius search without first requesting all nearest neighbors.
#'
#' The search uses unit vectors and Euclidean chord distances. Chord length
#' and angular distance are related by \eqn{d = 2\sin(\gamma/2)} and have
#' the same ordering for angles from 0 to 180 deg. Returned distances are
#' angular distances in radians, ordered from nearest to farthest.
#'
#' Sampling directions with missing angles are excluded. Query directions with
#' missing angles return empty neighborhoods. At polar angles of 0 or
#' `pi`, the second angle may be `NA`: the direction is still defined.
#' Returned indices always refer to the original sampling table.
#'
#' @param sampling_directions,query_directions data frames containing the angular
#'   columns selected by `angle_columns`, in radians, as returned by
#'   [add_directions()]. Other columns are ignored.
#' @param k positive integer-like numeric vector of length one, or `NULL`.
#'   Maximum number of neighbors. With `NULL`, retrieve all neighbors within
#'   `search_radius`, which must be supplied.
#' @param search_radius positive numeric vector of length one, at most 180,
#'   or `NULL`. Maximum angular distance in degrees. With `NULL`, only `k`
#'   limits the selection.
#' @param rule character vector of length one. With `"any"`, return the
#'   available neighbors. With `"all"`, return an empty neighborhood when
#'   fewer than `k` neighbors are available within `search_radius`. Ignored
#'   when `search_radius` or `k` is `NULL`.
#' @param angle_columns character vector of length two. Names of the polar
#'   and azimuthal angle columns, in that order. Defaults to
#'   `c("zenith", "azimuth")`. Polar angles must lie in `[0, pi]` and
#'   azimuthal angles in `[0, 2 * pi)`.
#' @inheritParams apply_locally
#'
#' @return list of length `nrow(query_directions)`, with one element per query:
#'   \describe{
#'     \item{`index`}{Integer vector of row positions in `sampling_directions`.}
#'     \item{`distance`}{Numeric vector of angular distances in radians.}
#'   }
#'   An empty neighborhood contains `integer(0)` and `numeric(0)`.
#'
#' @seealso [add_directions()], [calc_spherical_distance()], [apply_locally()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(10, 8),
#'   normalization_radius = 10,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(20, 16), geometry)
#' sampling_indices <- expand.grid(col = 1:20, row = 1:16)
#' m <- !is.na(g$off_axis) & !is.na(g$internal_azimuth)
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' query_indices <- data.frame(col = c(10, 14), row = c(8, 8))
#' sampling_directions <- add_directions(sampling_indices, g = g)
#' query_directions <- add_directions(query_indices, g = g)
#' neighbors <- nn_spherical(
#'   sampling_directions, query_directions, k = NULL, search_radius = 12,
#'   angle_columns = c("off_axis", "internal_azimuth"), parallel = FALSE
#' )
#' sampling_directions[neighbors[[1]]$index, ]
nn_spherical <- function(sampling_directions,
                         query_directions,
                         k = 1,
                         search_radius = NULL,
                         rule = "any",
                         parallel = TRUE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1,
                         angle_columns = c("zenith", "azimuth")) {
  .check_vector(k, "integerish", 1, sign = "positive", allow_null = TRUE)
  .check_vector(search_radius, "numeric", 1, sign = "positive",
                allow_null = TRUE)
  if (is.null(k) && is.null(search_radius)) {
    stop("Supply `search_radius` when `k = NULL`.", call. = FALSE)
  }
  if (!is.null(search_radius) && search_radius > 180) {
    stop("`search_radius` must not exceed 180 degrees.", call. = FALSE)
  }
  if (!is.null(search_radius) && !is.null(k)) {
    .assert_choice(rule, c("any", "all"))
  }
  .check_vector(angle_columns, "character", 2)
  if (anyDuplicated(angle_columns)) {
    stop("`angle_columns` must name two different columns.", call. = FALSE)
  }
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  sp <- .spherical_neighbor_vectors(sampling_directions, angle_columns)
  qp <- .spherical_neighbor_vectors(query_directions, angle_columns)
  sample_rows <- which(stats::complete.cases(sp))
  query_rows <- which(stats::complete.cases(qp))
  neighbors <- rep(list(list(index = integer(0), distance = numeric(0))),
                   nrow(query_directions))
  if (length(sample_rows) == 0 || length(query_rows) == 0) return(neighbors)
  sp <- sp[sample_rows, , drop = FALSE]
  qp <- qp[query_rows, , drop = FALSE]
  cores <- if (parallel) .cores(cores, logical, leave_free) else 1
  radius <- if (is.null(search_radius)) NULL else .degree2radian(search_radius)
  found <- nn_spherical_cpp(
    sp$x, sp$y, sp$z, qp$x, qp$y, qp$z,
    k = if (is.null(k)) 0L else min(k, nrow(sp)),
    search_radius = radius, cores = cores
  )
  for (i in seq_along(found)) {
    hit <- found[[i]]
    if (!is.null(search_radius) && !is.null(k) && rule == "all" &&
        length(hit$index) < k) next
    hit$index <- sample_rows[hit$index]
    neighbors[[query_rows[i]]] <- hit
  }
  neighbors
}

.spherical_neighbor_vectors <- function(directions, angle_columns) {
  .check_columns(directions, required_columns = angle_columns)
  polar <- directions[[angle_columns[1]]]
  azimuth <- directions[[angle_columns[2]]]
  .check_vector(polar, "numeric", allow_na = TRUE, name = angle_columns[1])
  .check_vector(azimuth, "numeric", allow_na = TRUE, name = angle_columns[2])
  if (any(polar < 0 | polar > pi, na.rm = TRUE) ||
      any(azimuth < 0 | azimuth >= 2 * pi, na.rm = TRUE)) {
    stop("Directions must be in [0, pi] and [0, 2 * pi) radians.",
         call. = FALSE)
  }
  x <- sin(polar) * cos(azimuth)
  y <- sin(polar) * sin(azimuth)
  z <- cos(polar)
  poles <- !is.na(polar) & polar %in% c(0, pi)
  x[poles] <- y[poles] <- 0
  data.frame(x = x, y = y, z = z)
}
