#' Remove statistical outliers in sky points
#'
#' @description
#' Remove sky points considered outliers relative to their local
#' neighbors in a user-specified variable.
#'
#' @param sampling_points `data.frame` with columns `"row"`, `"col"`, `"zen"`,
#'   `"azi"` (as in the output of [add_directions()]), and one additional
#'   numeric column with values to evaluate, identified by `col_id`.
#' @param col_id numeric or character vector of length one. The name or position
#'   of the column in `sampling_points` containing the values to evaluate.
#' @param search_radius numeric vector of length one. The maximum radius for
#'   searching k-nearest neighbors (KNN) in degrees.
#' @param detrend logical vector of length one. When `TRUE`, a detrending is
#'   attempted with [fit_cone_shaped_model()] with `method = zenith_n_azimuth`. If
#'   this attempt fails, the function reverts to `detrend = FALSE` and searches
#'   for outliers in the untransformed values. Detrending fails systematically
#'   when `k < 20`.
#'
#' @inheritParams calc_cv
#' @inheritParams extract_dn
#' @inheritParams add_directions
#' @inheritParams interpolate_spherical
#' @inheritParams is_outlier
#' @inheritParams apply_locally
#'
#' @details
#' Based on the Statistical Outlier Removal (SOR) filter from the
#' [PCL library](https://pointclouds.org/). Distances are computed on a spherical
#' surface. The number of neighbors is controlled by `k`, and `search_radius`
#' sets the maximum search radius (deg). If fewer than `k` neighbors are found
#' within that radius, the point is retained due to insufficient evidence for
#' removal. The decision criterion is from [is_outlier()].
#'
#' @return The retained points represented as a [data.frame], same as
#'   `sampling_points`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' gspec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), gspec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(gspec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#'
#' foo <- rem_outliers(sky_points, col_id = 5, k = 5, laxity = 2)
#'
#' plot(r)
#' points(sky_points$col, nrow(z) - sky_points$row, col = 2, pch = 20, cex = 1)
#' points(foo$col, nrow(z) - foo$row, col = 5, pch = 0)
#' }
rem_outliers <- function(sampling_points,
                         col_id,
                         k = 20,
                         search_radius = 20,
                         laxity = 2,
                         cutoff_side = "both",
                         detrend = FALSE,
                         parallel = FALSE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1) {

  col_id <- .resolve_required_column(col_id, sampling_points)
  .check_columns(sampling_points,
                 required_columns = c("col", "row", "zen", "azi", col_id),
                 reserved_columns = NULL)
  .check_vector(k, "integerish", 1)
  if (k < 3) stop("`k` must be at least 3.")
  .check_vector(detrend, "logical")

  .fn <- function(selected_points, ...) {
    .regular <- function(selected_points) {
      list(the_point_dn = selected_points[1, col_id],
           dns = selected_points[-1, col_id])
    }
    .detrend <- function(selected_points) {
      the_point <- selected_points[1, ]
      selected_points <- selected_points[-1, ]
      model <- fit_cone_shaped_model(selected_points, method = "zenith_n_azimuth")
      the_point_dn <- the_point[, col_id] - model$fun(the_point$zen, the_point$azi)
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

  neighbors <- nn_spherical(sampling_points, sampling_points,
                            k = k,
                            search_radius = search_radius,
                            rule = "all")
  i <- apply_locally(sampling_points,
                     sampling_points,
                     neighbors,
                     fun = .fn,
                     parallel = parallel,
                     cores = cores,
                     logical = logical,
                     leave_free = leave_free)
  i <- i$output
  i[is.na(i)] <- TRUE
  sampling_points[i, ]
}
