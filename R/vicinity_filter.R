#' Vicinity filter for planar and spherical coordinates
#'
#' Filters out nearby points based on a minimum distance threshold in either
#' planar (image space) or spherical (sky coordinates) systems.
#'
#' This function selects a subset of spatial points ensuring that no selected
#' points are closer than `min_dist`. Optionally, it can prioritize retention of
#' points with high values in a specific dimension such as image digital number.
#'
#' @inheritParams extract_rel_radiance
#' @param min_dist Numeric vector of length one. Minimum allowed distance
#'   between points. In degrees for spherical space or pixels for planar space.
#'   What space is used depend on the arguments `z` and `a`, if both are
#'   provided with the output of [zenith_image()] and [azimuth_image()],
#'   respectively, then the spherical space is used.
#' @param r [SpatRaster-class] or `NULL`. Set the dimension that will be
#'   prioritized when deciding to retein a sky point.
#' @param z [SpatRaster-class] built with [zenith_image()]  or `NULL`.
#' @param a [SpatRaster-class] built with [azimuth_image()]  or `NULL`.
#' @return The argument `sky_points` with fewer rows due to the removal of
#'   points closer each other than `min_dist`.
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
#' sky_points_p <- vicinity_filter(sky_points, r, min_dist = 100)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' points(sky_points_p$col, nrow(caim) - sky_points_p$row, col = 3, pch = 0)
#'
#' sky_points_s <- vicinity_filter(sky_points, r, z, a, min_dist = 30)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' points(sky_points_s$col, nrow(caim) - sky_points_s$row, col = 3, pch = 0)
#' }
vicinity_filter <- function(sky_points,
                            r = NULL,
                            z = NULL,
                            a = NULL,
                            min_dist = 3,
                            use_window = TRUE) {

  if (all(!is.null(r), !is.null(z), !is.null(z))) {
    .check_if_r_z_and_a_are_ok(r, z, a)
  }
  if (!is.null(r)) .is_single_layer_raster(r, "r")
  if (any(!is.null(z), !is.null(a))) {
    .is_single_layer_raster(z, "z")
    .is_single_layer_raster(a, "a")
    terra::compareGeom(z, a)
  }
  stopifnot(is.numeric(min_dist))
  stopifnot(length(min_dist) == 1)
  stopifnot(is.data.frame(sky_points))
  stopifnot(ncol(sky_points) == 2)
  stopifnot(nrow(sky_points) > 1)

  .dist_spherical <- function(zenith_azimuth) {
    n <- nrow(zenith_azimuth)
    m <- matrix(0, n, n)

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        d <- calc_spherical_distance (zenith_azimuth[i, 1],
                                      zenith_azimuth[i, 2],
                                      zenith_azimuth[j, 1],
                                      zenith_azimuth[j, 2])
        m[i, j] <- d
        m[j, i] <- d
      }
    }
    dimnames(m) <- list(rownames(zenith_azimuth), rownames(zenith_azimuth))
    m
  }


  if (!is.null(r)) {
    ord <- extract_dn(r, sky_points, use_window = use_window)[,3]
    ord <- order(ord, decreasing = TRUE)
    sky_points <- sky_points[ord, ]
  }

  if (is.null(z)) {
    coords <- sky_points[, c("col", "row")]
    dists <- as.matrix(stats::dist(coords))
  } else {
    sky_points <- extract_dn(c(z, a), sky_points, use_window = FALSE)
    coords <- sky_points[, c(3, 4)] %>% .degree2radian()
    dists <- .dist_spherical(coords)
    min_dist <- .degree2radian(min_dist)
  }

  n <- nrow(coords)
  selected <- logical(n)
  considered <- logical(n)
  to_consider <- seq_len(n)
  while (length(to_consider) > 0) {
    i <- to_consider[1]
    selected[i] <- TRUE
    considered[i] <- TRUE
    too_close <- which(dists[i, ] <= min_dist)
    considered[too_close] <- TRUE
    to_consider <- which(!considered)
  }
  sky_points[selected, c("row", "col"), drop = FALSE]
}
