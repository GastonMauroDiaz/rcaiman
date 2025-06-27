#' Statistical outlier removal (SOR) filter
#'
#' Statistical outlier removal (SOR) filter
#'
#' This algorithm is based on the homonymous filter from the [PCL
#' library](https://pointclouds.org/). Distances are computed on a spherical
#' surface and expressed in degrees to avoid distortions due to projection. The
#' number of neighbors used for evaluation is controlled by the `k` argument,
#' while the `chi_max` argument sets the maximum search radius for finding these
#' neighbors. Points are projected onto a unit-radius sphere, similar to the use
#' of relative radius in image mapping. The spherical distance is then
#' calculated, and points farther than chi_max are excluded from the neighbor set.
#' If an insufficient number of neighbors are found within chi_max, the point is
#' retained due to a lack of evidence for removal. The decision criterion
#' follows \insertCite{Leys2013;textual}{rcaiman}:
#'
#' \eqn{M - laxity \times MAD < x_i < M + laxity \times MAD}
#'
#' where \eqn{x_i} is the value associated with a given sky point, \eqn{M} and
#' \eqn{MAD} are the median and median absolute deviation, respectively,
#' computed from the values associated with the neighbors of \eqn{x_i}, and
#' \eqn{laxity} is the user-defined threshold.
#'
#' The argument `cutoff_side` controls which side of the inequality is
#' tested—either one side or both.
#'
#' @inheritParams extract_dn
#' @inheritParams sky_grid_segmentation
#' @inheritParams interpolate_planar
#' @param chi_max Numeric vector of length one. The maximum radius for searching
#'   k-nearest neighbors (knn). Points are projected onto a unit-radius sphere,
#'   similar to the use of relative radius in image mapping. The spherical
#'   distance is then calculated and used to filter out points farther than
#'   `chi_max`.The distance is expressed in degrees. If an insufficient number of
#'   neighbors are found within the search radius, the point is retained due to
#'   a lack of evidence for removal.
#' @param laxity Numeric vector of length one. See
#'   \insertCite{Leys2013;textual}{rcaiman}.
#' @param cutoff_side Character vector of length one. See
#'   \insertCite{Leys2013;textual}{rcaiman}.
#' @param trend Numeric vector of length one or `NULL`. Specifies the order of
#'   the surface fitted to the set of points found by the KNN algorithm in order
#'   to remove trends such as smooth gradients. The maximum allowed order is
#'   three. If `NULL` (default), no attempt is made to remove a trend.
#'
#' @return The argument `sky_points` with fewer rows due to outlier removal.
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
#' sky_points <- sor_filter(sky_points, r, z, a, k = 10, chi_max = 20, laxity = 2,
#'                 cutoff_side = "left")
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 3, pch = 0)
#'
#' }
sor_filter <- function(sky_points, r, z, a,
                       k = 20,
                       chi_max = 20,
                       laxity = 2,
                       cutoff_side = "both",
                       use_window = TRUE,
                       trend = NULL) {

  .check_if_r_z_and_a_are_ok(r, z, a)
  stopifnot(length(k) == 1)
  stopifnot(.is_whole(k))
  stopifnot(k >= 3)
  stopifnot(length(chi_max) == 1)
  stopifnot(is.numeric(chi_max))
  chi_max <- .degree2radian(chi_max)

  sky_points <- extract_dn(c(z, a, r), sky_points, use_window = TRUE)
  names(sky_points)[3:5] <- c("z", "a", "dn")

  .calculate_sor <- function(i) {
    sky_points[, c("z", "a")] <- .degree2radian(sky_points[, c("z", "a")])
    spherical_distance <- calc_spherical_distance (sky_points$z,
                                                   sky_points$a,
                                                   sky_points[i, "z"],
                                                   sky_points[i, "a"])
    order_idx <- order(spherical_distance)
    sorted_distance <- spherical_distance[order_idx][2:(k + 1)]
    if (!is.null(trend)) stopifnot(trend <= 6)
    tryCatch(
      if (all(sorted_distance <= chi_max)) {
        if (!is.null(trend)) {
          the_point <- sky_points[order_idx[1], c("row", "col", "dn")]
          sky_points <- sky_points[order_idx[2:(k + 1)], c("row", "col", "dn")]
          xy <- terra::cellFromRowCol(r, sky_points$row, sky_points$col) %>%
            terra::xyFromCell(r, .)
          fit <- spatial::surf.ls(x = xy[, 1],
                                  y = xy[, 2],
                                  z = sky_points[, "dn"],
                                  np = trend)
          pred <- predict(fit, xy[,1], xy[,2])
          dns <- fit$z - pred
          xy <- terra::cellFromRowCol(r, the_point$row, the_point$col) %>%
            terra::xyFromCell(r, .)
          the_point_dn <- the_point[, "dn"] - predict(fit, xy[,1], xy[,2])
        } else {
          the_point_dn <- sky_points[order_idx[1], "dn"]
          dns <- sky_points[order_idx[2:(k + 1)], "dn"]
        }
        return(c(the_point_dn,
                 stats::median(dns, na.rm = TRUE),
                 stats::mad(dns, na.rm = TRUE)))
      } else {
        return(c(NA, NA, NA))
      },
      error = function(e) stop("Try other order of magnitud for 'trend' or look
                               for sky points below the horizon or on it.")
    )
  }

  result <- lapply(seq_len(nrow(sky_points)), .calculate_sor) %>% unlist() %>%
      matrix(., ncol = 3, byrow = TRUE)

  value <- result[, 1]
  central_tendency <- result[, 2]
  dispersion <- result[, 3]

  deviation <- (value - central_tendency) / dispersion
  i <- switch(cutoff_side,
              right = deviation <= laxity,
              left = deviation >= -laxity,
              both = abs(deviation) <= laxity)

  i[is.na(i)] <- TRUE
  sky_points[i, c("row", "col") ]
}
