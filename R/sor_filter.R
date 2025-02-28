#' Statistical outlier removal filter
#'
#' @inheritParams extract_rl
#' @param r [SpatRaster-class]. An image with the same raster grid as the one
#'   from which `sky_points` was obtained.
#' @inheritParams interpolate_sky_points
#' @param rmax Numeric vector of length one. The maximum radius for searching
#'   k-nearest neighbors (knn). Points are projected onto a unit-radius sphere,
#'   similar to the use of relative radius in image mapping. The spherical
#'   distance is then calculated and used to filter out points farther than
#'   `rmax`. The distance is expressed in degrees.
#' @param thr Numeric vector of length one. See
#'   \insertCite{Leys2013;textual}{rcaiman}.
#' @param cutoff_side Character vector of length one. See
#'   \insertCite{Leys2013;textual}{rcaiman}.
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
#' bin <- bin & mask_hs(z, 0, 80)
#' g <- sky_grid_segmentation(z, a, 5, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g,
#'                                  dist_to_plant = 3,
#'                                  min_raster_dist = 10)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' sky_points <- extract_rl(r, z, a, sky_points, no_of_points = NULL)
#' i <- sor_filter(sky_points$sky_points, r, rmax = 10, thr = 1,
#'                 cutoff_side = "left")
#' sky_points <- sky_points$sky_points[i, c("row", "col")]
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' }
sor_filter <- function(sky_points, r,
                       k = 20,
                       rmax = 20,
                       thr = 2,
                       cutoff_side = "both") {
  stopifnot(k >= 3)
  stopifnot(length(rmax) == 1)
  stopifnot(is.numeric(rmax))
  rmax <- .degree2radian(rmax)

  ds <- extract_dn(r, sky_points[, c("row", "col")])

  calculate_sor <- function(i) {
    spherical_distance <- .calc_spherical_distance(sky_points$z,
                                                   sky_points$a,
                                                   sky_points[i, "z"],
                                                   sky_points[i, "a"],
                                                   radians = FALSE)
    order_idx <- order(spherical_distance)
    sorted_distance <- spherical_distance[order_idx][2:(k + 1)]
    m <- sorted_distance <= rmax

    if (sum(m) >= 3) {
      u <- ds[order_idx[2:(k + 1)], 3]
      return(c(stats::median(u[m], na.rm = TRUE), stats::mad(u[m], na.rm = TRUE)))
    } else {
      return(c(NA, NA))
    }
  }

  result <- lapply(seq_len(nrow(sky_points)), calculate_sor) %>% unlist() %>%
      matrix(., ncol = 2, byrow = TRUE)

  central_tendency <- result[, 1]
  dispersion <- result[, 2]

  deviation <- (ds[, 3] - central_tendency) / dispersion
  i <- switch(cutoff_side,
              right = deviation < thr,
              left = deviation > -thr,
              both = abs(deviation) < thr)

  i[is.na(i)] <- TRUE
  i
}
