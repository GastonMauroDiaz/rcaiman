#' Statistical outlier removal (SOR) filter
#'
#' Statistical outlier removal (SOR) filter
#'
#' This algorithm is based on the homonymous filter from the [PCL
#' library](https://pointclouds.org/). Distances are computed on a spherical
#' surface and expressed in degrees to avoid distortions due to projection. The
#' number of neighbors used for evaluation is controlled by the `k` argument, while
#' the `rmax` argument sets the maximum search radius for finding these neighbors.
#' Points are projected onto a unit-radius sphere, similar to the use of relative
#' radius in image mapping. The spherical distance is then calculated, and points
#' farther than rmax are excluded from the neighbor set. If an insufficient
#' number of neighbors are found within rmax, the point is retained due to a
#' lack of evidence for removal. The decision criterion follows
#' \insertCite{Leys2013;textual}{rcaiman}:
#'
#' \eqn{M - thr \times MAD < x_i < M + thr \times MAD}
#'
#' where \eqn{x_i} is the value associated with a given sky point, \eqn{M} and
#' \eqn{MAD} are the median and median absolute deviation, respectively,
#' computed from the values associated with the neighbors of \eqn{x_i}, and
#' \eqn{thr} is the user-defined threshold.
#'
#' The argument `cutoff_side` controls which side of the inequality is
#' tested—either one side or both.
#'
#' @inheritParams extract_dn
#' @inheritParams ootb_mblt
#' @inheritParams interpolate_sky_points
#' @param rmax Numeric vector of length one. The maximum radius for searching
#'   k-nearest neighbors (knn). Points are projected onto a unit-radius sphere,
#'   similar to the use of relative radius in image mapping. The spherical
#'   distance is then calculated and used to filter out points farther than
#'   `rmax`.The distance is expressed in degrees. If an insufficient number of
#'   neighbors are found within the search radius, the point is retained due to
#'   a lack of evidence for removal.
#' @param thr Numeric vector of length one. See
#'   \insertCite{Leys2013;textual}{rcaiman}.
#' @param cutoff_side Character vector of length one. See
#'   \insertCite{Leys2013;textual}{rcaiman}.
#'
#' @return The argument `sky_points` with fewer rows due to outlier removal.
#'
#' @family Tool Functions
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' bin <- bin & select_sky_vault_region(z, 0, 80)
#' g <- sky_grid_segmentation(z, a, 5, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g,
#'                                  dist_to_black = 3)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- sor_filter(sky_points, r, z, a, k = 10, rmax = 20, thr = 2,
#'                 cutoff_side = "left")
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 3, pch = 0)
#'
#' }
sor_filter <- function(sky_points, r, z, a,
                       k = 20,
                       rmax = 20,
                       thr = 2,
                       cutoff_side = "both",
                       use_window = TRUE) {

  .check_if_r_z_and_a_are_ok(r, z, a)
  stopifnot(length(k) == 1)
  stopifnot(.is_whole(k))
  stopifnot(k >= 3)
  stopifnot(length(rmax) == 1)
  stopifnot(is.numeric(rmax))
  rmax <- .degree2radian(rmax)

  sky_points <- extract_dn(c(z, a, r), sky_points, use_window = TRUE)
  names(sky_points)[3:5] <- c("z", "a", "dn")

  .calculate_sor <- function(i) {
    spherical_distance <- .calc_spherical_distance(
                                      sky_points$z %>% .degree2radian(),
                                      sky_points$a %>% .degree2radian(),
                                      sky_points[i, "z"] %>% .degree2radian(),
                                      sky_points[i, "a"] %>% .degree2radian(),
                                      radians = TRUE)
    order_idx <- order(spherical_distance)
    sorted_distance <- spherical_distance[order_idx][2:(k + 1)]
    tryCatch(
    if (all(sorted_distance <= rmax)) {
      dns <- sky_points[order_idx[2:(k + 1)], "dn"]
      return(c(stats::median(dns, na.rm = TRUE), stats::mad(dns, na.rm = TRUE)))
    } else {
      return(c(NA, NA))
    },
      error = function(e) stop("There are sky points below the horizon or on it.")
    )
  }

  result <- lapply(seq_len(nrow(sky_points)), .calculate_sor) %>% unlist() %>%
      matrix(., ncol = 2, byrow = TRUE)

  central_tendency <- result[, 1]
  dispersion <- result[, 2]

  deviation <- (sky_points[, "dn"] - central_tendency) / dispersion
  i <- switch(cutoff_side,
              right = deviation <= thr,
              left = deviation >= -thr,
              both = abs(deviation) <= thr)

  i[is.na(i)] <- TRUE
  sky_points[i, c("row", "col") ]
}
