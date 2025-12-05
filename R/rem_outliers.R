#' Remove statistical outliers in sky points
#'
#' @description
#' Remove sky points considered outliers relative to their local
#' neighbors in a user-specified variable.
#'
#' @param angular_radius numeric vector of length one. The maximum radius for
#'   searching k-nearest neighbors (KNN) in degrees.
#' @param trend numeric vector of length one or `NULL`. Zero to three. Specifies
#'   the order of the polynomial surface fitted to the neighbors to account for
#'   spatial trends. Use NULL (default) to skip detrending.
#'
#' @inheritParams extract_dn
#' @inheritParams skygrid_segmentation
#' @inheritParams interpolate_spherical
#' @inheritParams is_outlier
#' @inheritParams apply_by_direction
#'
#' @details
#' Based on the Statistical Outlier Removal (SOR) filter from the
#' [PCL library](https://pointclouds.org/). Distances are computed on a spherical
#' surface. The number of neighbors is controlled by `k`, and `angular_radius`
#' sets the maximum search radius (deg). If fewer than `k` neighbors are found
#' within that radius, the point is retained due to insufficient evidence for
#' removal. The decision criterion is from [is_outlier()].
#'
#' @note This function assumes that `sky_points` and the
#' [terra::SpatRaster-class] objects refer to the same image geometry. No checks
#' are performed.
#'
#' @return The retained points represented as a [data.frame] with columns `row`
#'   and `col`, same as `sky_points`.
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
#' bin <- binarize_by_region(r, ring_segmentation(z, 30),
#'                           method = "thr_isodata")
#' bin <- bin & select_sky_region(z, 0, 80)
#' g <- skygrid_segmentation(z, a, 5, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g,
#'                                  dist_to_black = 3)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- rem_outliers(sky_points, r, z, a,
#'                                  k = 5,
#'                                  angular_radius = 20,
#'                                  laxity = 2,
#'                                  cutoff_side = "left")
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 3, pch = 0)
#' }
rem_outliers <- function(sky_points, r, z, a,
                         k = 20,
                         angular_radius = 20,
                         laxity = 2,
                         cutoff_side = "both",
                         use_window = TRUE,
                         trend = NULL,
                         parallel = FALSE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1) {

  .check_sky_points(sky_points)
  .check_r_z_a_m(r, z, a, r_type = "single")
  if (k < 3) stop("`k` must be at least 3.")
  .check_vector(use_window, "logical", 1)
  .check_vector(trend, "integerish", 1, allow_null = TRUE,  sign = "positive")
  if (!is.null(trend) && trend > 3) stop("`trend` must be 1-3 or NULL.")

  sky_points <- extract_dn(c(z, a, r), sky_points, use_window = use_window)
  names(sky_points)[3:5] <- c("z", "a", "dn")

  .fn <- function(sky_points) {
    .regular <- function(sky_points) {
      list(the_point_dn = sky_points[1, "dn"],
           dns = sky_points[-1, "dn"])
    }
    .detrend <- function(sky_points) {
      the_point <- sky_points[1, ]
      sky_points <- sky_points[-1, ]
      xy <- terra::cellFromRowCol(r, sky_points$row, sky_points$col) %>%
        terra::xyFromCell(r, .)
      fit <- tryCatch(spatial::surf.ls(x = xy[, 1],
                              y = xy[, 2],
                              z = sky_points[, "dn"],
                              np = trend),
                      error = function(e) {
                        warning("Detrending failed")
                        stop()
                      })
      pred <- predict(fit, xy[,1], xy[,2])
      dns <- fit$z - pred
      xy <- terra::cellFromRowCol(r, the_point$row, the_point$col) %>%
        terra::xyFromCell(r, .)
      the_point_dn <- the_point[, "dn"] - predict(fit, xy[,1], xy[,2])
      list(the_point_dn = the_point_dn, dns = dns)
    }
    l <- if (!is.null(trend)) {
      tryCatch(.detrend(sky_points), error = function(e) .regular(sky_points))
    } else {
      .regular(sky_points)
    }
    !is_outlier(l$the_point_dn, l$dns, laxity = laxity, cutoff_side = cutoff_side)
  }

  .fn(sky_points)

  i <- apply_local_spherical(sky_points, NULL, z, a, k, angular_radius,
                             rule = "all", fun = .fn,
                             parallel = parallel,
                             cores = cores,
                             logical = logical,
                             leave_free = leave_free)
  i <- i$output
  i[is.na(i)] <- TRUE
  sky_points[i, c("row", "col") ]
}
