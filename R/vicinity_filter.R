#' Vicinity filter for planar and spherical coordinates
#'
#' Filters out nearby points based on a minimum distance threshold in either
#' planar (image space) or spherical (sky coordinates) systems.
#'
#' This function selects a subset of spatial points ensuring that no selected
#' points are closer than `thr`. Useful for obtaining well-distributed samples.
#'
#' @inheritParams extract_rel_radiance
#' @param thr Numeric vector of length one. Minimum allowed distance between
#'   points. In degrees for spherical space or pixels for planar space.
#' @param space Numeric vector of length one. Either `"planar"` or
#'   `"spherical"`.
#' @param prefer Character vector of length one or `NULL`. Optional. Name of
#'   column whose higher values are to be prioritized during selection.
#'
#' @return The argument `sky_points` with fewer rows due to the removal of
#'   points closer each other than `thr`.
#'
#' @family Tool Functions
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
#'
#' sky_points_p <- vicinity_filter(sky_points, "planar", 100)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' points(sky_points_p$col, nrow(caim) - sky_points_p$row, col = 3, pch = 0)
#'
#' rr <- extract_rel_radiance(r, z, a, sky_points)
#' sky_points_s <- vicinity_filter(rr$sky_points, "spherical", 30)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' points(sky_points_s$col, nrow(caim) - sky_points_s$row, col = 3, pch = 0)
#' }
vicinity_filter <- function(sky_points, space, thr, prefer = NULL) {
  stopifnot(is.numeric(thr))
  stopifnot(length(thr) == 1)
  stopifnot(space %in% c("planar", "spherical"))
  stopifnot(nrow(sky_points) > 1)

  .dist_spherical <- function(zenith_azimuth) {
    n <- nrow(zenith_azimuth)
    m <- matrix(0, n, n)

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        d <- .calc_spherical_distance(zenith_azimuth[i, 1],
                                      zenith_azimuth[i, 2],
                                      zenith_azimuth[j, 1],
                                      zenith_azimuth[j, 2],
                                      radians = TRUE)
        m[i, j] <- d
        m[j, i] <- d
      }
    }
    dimnames(m) <- list(rownames(zenith_azimuth), rownames(zenith_azimuth))
    m
  }

  if (!is.null(prefer)) {
    stopifnot(prefer %in% colnames(sky_points))
    ord <- order(sky_points[[prefer]], decreasing = TRUE)
    sky_points <- sky_points[ord, , drop = FALSE]
  }

  if (space == "planar") {
    coords <- sky_points[, c("col", "row")]
    dists <- as.matrix(stats::dist(coords))
  } else {
    coords <- sky_points[, c("z", "a")] %>% .degree2radian()
    dists <- .dist_spherical(coords)
    thr <- .degree2radian(thr)
  }

  n <- nrow(coords)
  selected <- logical(n)
  considered <- logical(n)
  to_consider <- seq_len(n)
  while (length(to_consider) > 0) {
    i <- to_consider[1]
    selected[i] <- TRUE
    considered[i] <- TRUE
    too_close <- which(dists[i, ] <= thr)
    considered[too_close] <- TRUE
    to_consider <- which(!considered)
  }
  sky_points[selected, , drop = FALSE]
}
