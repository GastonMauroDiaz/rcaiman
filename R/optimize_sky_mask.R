#' Optimize sky mask selection
#'
#' @description
#' Select the binarized image that provides the best trade-off between
#' radiometric consistency of sky pixels and spatial coverage.
#'
#' @details
#' Each element in `bin_list` is evaluated as a candidate sky mask. For each
#' mask, sky points are generated using [fibonacci_points()] at a fixed
#' angular resolution of 1 deg, ensuring high density sampling.
#'
#' Two criteria are computed:
#' \itemize{
#'   \item *radiometric consistency*, estimated from the local variability provided by `cv`;
#'   \item *coverage*, the fraction of equal-area cells containing at least one
#'     sampled sky point.
#' }
#'
#' Radiometric consistency final value for each mask is the median of CV values
#' across sampled points. Lower values indicate more homogeneous (and thus more
#' reliable) sky samples. This assumes that true sky regions exhibit
#' smoother radiance patterns than canopy elements. Uniform non-sky regions
#' (e.g., deep shadows) may also produce low variability, but such cases are
#' typically associated with poor binarizations and are not expected to be even
#' considered in the candidate set.
#'
#' Coverage is computed with [assess_sampling_uniformity()].
#'
#' Both components are scaled to \eqn{[0,1]} using [normalize_min_max()].
#'
#' After normalization, the Coverage-Consistency Metric is defined as:
#'
#' \deqn{
#' \mathrm{CCM} =
#' (1 - \mathrm{coverage}_n) \cdot w +
#' \mathrm{consistency}_n \cdot (1 - w)
#' }
#'
#' where \eqn{n} is the nth binary mask being tested. Lower values indicate
#' better trade-offs between high coverage and low radiometric variability.
#' Restrictive masks may improve consistency but reduce coverage, whereas
#' permissive masks increase coverage at the cost of including mixed or
#' non-sky pixels. The maximum achievable coverage (before normalization)
#' depends on canopy openness (see [compute_canopy_openness()]), as sampling
#' becomes more efficient under less obstructed sky conditions.
#'
#' @param cv [terra::SpatRaster-class] generated with [cv_image()].
#' @param bin_list list of logical [terra::SpatRaster-class] objects of one
#'   layer.
#' @param w numeric vector of length one. Weight controlling the balance between
#'   coverage and consistency (see *Details*).
#'
#' @inheritParams add_directions
#' @inheritParams assess_sampling_uniformity
#'
#' @return Logical [terra::SpatRaster-class]. Optimal sky mask selected from `bin_list`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' thr <- apply_by_direction(r, z, a, m)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' seg <- equal_area_segmentation(z, a, 1000)
#' sky_points <- sample_sky_points(r, bin, seg, dist_to_black = 1)
#' sky_points <- add_directions(sky_points, z, a)
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- rem_outliers(
#'   sky_points,
#'   col_id = "Cyan",
#'   k = 20,
#'   search_radius = 20,
#'   laxity = 2,
#'   cutoff_side = "left"
#' )
#' cv <- cv_image(r)
#' sky_points <- extract_dn(cv, sky_points, use_window = FALSE)
#' i <- is_outlier(sky_points$cv, sky_points$cv, laxity = 1, cutoff_side = "right")
#' sky_points <- sky_points[!i,]
#'
#' sky <- sky_via_local_models(z, a, sky_points, col_id = "Cyan")
#' plot(sky)
#' plot(r/sky)
#'
#' bin_list <- lapply(list(0.5, 0.7, 0.9), function(w) {
#'   bin <- binarize_with_thr(r/sky, w)
#'   bin <- rem_isolated_black_pixels(bin)
#'   lapply(1:9, function(x) {
#'     bin <- grow_black(bin, x)
#'     names(bin) <- paste0("w_", w, "_grow_", x)
#'     bin
#'   })
#' })
#' bin_list <- (unlist(bin_list, recursive = FALSE))
#'
#' equal_area_seg <- equal_area_segmentation(z, a, 200)
#'
#' bin <- optimize_sky_mask(cv, z, a, equal_area_seg, bin_list, w = 0.5)
#' bin
#' plot(bin)
#' }
optimize_sky_mask <- function(cv, z, a,
                              equal_area_seg,
                              bin_list,
                              w = 0.5) {

  .assert_single_layer(cv, squared = TRUE)
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .assert_same_geom(z, cv)
  .assert_equal_area_seg(equal_area_seg)

  lapply(seq_along(bin_list),
         function(i) {
           .assert_logical_mask(bin_list[[i]])
           .assert_same_geom(bin_list[[i]], cv)
           }
         )

  .check_vector(w, "numeric", 1, sign = "nonnegative")
  if (w > 1) stop("`w` cannot be greater than one.")

  bin_list_index <- seq_along(bin_list)

  metrics <- lapply(seq_along(bin_list_index), function(i) {
    idx <- bin_list_index[i]
    sky_points <- fibonacci_points(z, a, 1, m = bin_list[[idx]])
    sky_points <- extract_dn(cv, sky_points, use_window = FALSE)

    acc <- stats::median(sky_points$cv)
    cov <- assess_sampling_uniformity(
      sky_points,
      equal_area_seg
    )$coverage

    list(acc = acc,
         cov = cov,
         bin_list_index = idx)
  })

  .fun <- function(id = "acc") {
      vapply(metrics, function(x) x[[id]], numeric(1))
  }

  acc_vec <- .fun("acc")
  cov_vec <- .fun("cov")

  acc_n <- normalize_min_max(acc_vec)
  cov_n <- normalize_min_max(cov_vec)

  cam_vec <- (1 - cov_n) * w + acc_n * (1 - w)
  best <- which.min(cam_vec)

  bin_list_index_vec <- .fun("bin_list_index")
  i <- bin_list_index_vec[best]

  bin_list[[i]]
}
