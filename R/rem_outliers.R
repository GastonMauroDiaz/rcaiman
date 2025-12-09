#' Remove statistical outliers in sky points
#'
#' @description
#' Remove sky points considered outliers relative to their local
#' neighbors in a user-specified variable.
#'
#' @param angular_radius numeric vector of length one. The maximum radius for
#'   searching k-nearest neighbors (KNN) in degrees.
#' @param detrend logical vector of length one. When `TRUE`, a detrending is
#'   attempted with [fit_coneshaped_model()] with `method = zenith_n_azimuth`. If
#'   this attempt fails, the function reverts to `detrend = FALSE` and searches
#'   for outliers in the untransformed values. Detrending fails systematically
#'   when `k < 20`.
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
#' seg <- equalarea_segmentation(z, a, 1000)
#' sky_points <- sample_sky_points(r, bin, seg,
#'                                 dist_to_black = 3)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points2 <- rem_outliers(sky_points, r, z, a,
#'                             k = 20,
#'                             angular_radius = 30,
#'                             laxity = 2,
#'                             cutoff_side = "left",
#'                             detrend = TRUE)
#' points(sky_points2$col, nrow(caim) - sky_points2$row, col = 3, pch = 0)
#' }
rem_outliers <- function(sky_points, r, z, a,
                         k = 20,
                         angular_radius = 20,
                         laxity = 2,
                         cutoff_side = "both",
                         use_window = TRUE,
                         detrend = FALSE,
                         parallel = FALSE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1) {

  .check_sky_points(sky_points)
  .check_r_z_a_m(r, z, a, r_type = "single")
  if (k < 3) stop("`k` must be at least 3.")
  .check_vector(use_window, "logical", 1)
  .check_vector(detrend, "logical")

  sky_points <- extract_dn(c(z, a, r), sky_points, use_window = use_window)
  names(sky_points)[3:5] <- c("z", "a", "rr")

  .fn <- function(selected_points) {
    .regular <- function(selected_points) {
      list(the_point_dn = selected_points[1, "rr"],
           dns = selected_points[-1, "rr"])
    }
    .detrend <- function(selected_points) {
      the_point <- selected_points[1, ]
      selected_points <- selected_points[-1, ]
      model <- fit_coneshaped_model(selected_points, method = "zenith_n_azimuth")
      the_point_dn <- the_point[, "rr"] - model$fun(the_point$z, the_point$a)
      dns <- stats::residuals(model$model)
      list(the_point_dn = the_point_dn, dns = dns)
    }
    l <- if (detrend) {
      tryCatch(.detrend(selected_points), error = function(e) .regular(selected_points))
    } else {
      .regular(selected_points)
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
