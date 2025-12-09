#' Apply a function locally in the spherical domain
#'
#' @description
#' Evaluate a user-supplied function locally using neighboring sky points
#' defined on the hemispherical domain.
#'
#' @details
#' The `query_points` argument lets the user evaluate `fun` at a custom set of
#' locations. For each point in `query_points`, neighbors are searched within
#' `sky_points`. When `query_points` is `NULL` (default), the evaluation is
#' performed for every point in `sky_points`. Neighbor selection is based on
#' spherical distance, see [calc_spherical_distance()]. The function `fun` is
#' then applied to the corresponding subset of neighboring sky points.
#'
#' The input passed to `fun` is a **subset of rows of `sky_points`**, therefore,
#' `fun` must extract and process the column(s) of interest, e.g.: `function(df)
#' mean(df$dn)`.
#'
#' The angular neighborhood is defined by two parameters:
#' \itemize{
#'   \item `k`: number of nearest neighbors used by `fun`.
#'   \item `angular_radius`: maximum allowed spherical distance (deg).
#' }
#'
#' Argument `rule` controls how neighbor availability is interpreted:
#' \itemize{
#'   \item `"any"`: all neighbors found within `angular_radius` are passed to `fun`
#'                 (even if fewer than `k`).
#'   \item `"all"`: if fewer than `k` neighbors lie within the radius, the query
#'                 point is assigned `NA`.
#' }
#'
#' @param sky_points data.frame with columns `row` and `col`. Source points used
#'   as candidate neighbors. See [sample_sky_points()] and [extract_dn()].
#' @param query_points optional data.frame with columns `row` and
#'   `col` with the points at which `fun` should be evaluated.
#' @param fun function. Must accept a data.frame (subset of `sky_points`) and
#'   return a numeric scalar.
#' @param rule character vector of length one. `"any"` or `"all"`. See *Details*.
#'
#' @inheritParams skygrid_segmentation
#' @inheritParams interpolate_spherical
#' @inheritParams apply_by_direction
#'
#' @seealso
#' [skygrid_centers()] and [fibonacci_points()] for producing `query_points`.
#'
#' @return
#' A data.frame equal to `query_points` with one additional column, named
#' `"out"`, containing the result of applying `fun` at each query point.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' sky_points <- fibonacci_points(z, a, 3)
#' display_caim(caim, sky_points = sky_points)
#' sky_points <- extract_dn(r, sky_points, use_window = FALSE)
#' head(sky_points)
#' eval_points <- apply_local_spherical(
#'   sky_points,
#'   fibonacci_points(z, a, 10),
#'   z,
#'   a,
#'   k = 500,
#'   angular_radius = 30,
#'   rule = "any",
#'   fun = function(df) thr_isodata(df$Blue)
#' )
#'
#' thr <- triangulate(eval_points[!is.na(eval_points$out),], r, col_id = 3)
#' plot(thr)
#' plot(binarize_with_thr(r, thr))
#' }
apply_local_spherical <- function(sky_points,
                                  query_points = NULL,
                                  z,
                                  a,
                                  k = 20,
                                  angular_radius = 20,
                                  rule = "any",
                                  fun,
                                  parallel = TRUE,
                                  cores = NULL,
                                  logical = TRUE,
                                  leave_free = 1
) {

  if (!is.null(query_points)) .check_sky_points(query_points)
  required_cols <- c("row", "col")
  if (!all(required_cols %in% names(sky_points))) {
    stop(sprintf("`sky_points` must at least contain columns %s.",
                 paste(sprintf('"%s"', required_cols), collapse = ", ")))
  }
  .check_r_z_a_m(r = NULL, z, a)
  .check_vector(k, "integerish", 1, sign = "positive")
  .check_vector(angular_radius, "numeric", 1, sign = "positive")
  .assert_choice(rule, c("any", "all"))
  if (!is.function(fun))
    stop("`fun` must be a function returning a numeric vector of length one.")
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }

  # convert radius to radians
  angular_radius <- .degree2radian(angular_radius)

  .extract_angles_radians <- function(x) {
    x <- extract_dn(c(z, a), x[c("row", "col")], use_window = FALSE)
    x <- .degree2radian(x[c("Zenith image", "Azimuth image")])
    names(x) <- c("z", "a")
    x
  }

  s_angles <- .extract_angles_radians(sky_points)

  if (!is.null(query_points)) {
    q_angles <- .extract_angles_radians(query_points)
    n <- nrow(query_points)
    process_one <- function(i) {
      d <- calc_spherical_distance(s_angles$z, s_angles$a,
                                   q_angles[i, "z"],
                                   q_angles[i, "a"])

      inside <- which(d <= angular_radius)

      if (rule == "all" && length(inside) < k) return(NA)
      if (length(inside) < 1) return(NA)

      nn <- inside[order(d[inside])[1:k]]
      nn <- nn[!is.na(nn)]

      fun(sky_points[nn,])
    }
  } else {
    n <- nrow(sky_points)
    process_one <- function(i) {
      d <- calc_spherical_distance(s_angles$z, s_angles$a,
                                   s_angles[i, "z"],
                                   s_angles[i, "a"])

      # neighbors sorted by distance
      order_idx <- order(d)

      # get k neighbors excluding self
      neighbor_idx <- order_idx[2:(k + 1)]
      neighbor_dist <- d[neighbor_idx]

      # select neighbors within radius
      nn <- neighbor_dist <= angular_radius
      if (rule == "all" && sum(nn) < k) return(NA)
      if (sum(nn) < 1) return(NA)
      nn <- nn[!is.na(nn)]

      fun(sky_points[c(order_idx[1], neighbor_idx[nn]), ])
    }
  }

  output <- if (parallel) {
    .with_cluster(cores, {
      # Only to avoid note from check, code is OK without this line.
      i <- NA

      foreach::foreach(i = seq_len(n), .combine = c) %dopar% {
        list(process_one(i))
      }
    })
  } else {
    lapply(seq_len(n), process_one)
  }
  output <- unlist(output)
  if (!is.null(query_points)) {
    return(cbind(query_points, output))
  } else {
    return(cbind(sky_points, output))
  }
}
