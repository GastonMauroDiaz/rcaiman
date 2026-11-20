#' Optimize sky mask selection
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Select a candidate sky mask by balancing local radiometric variability and
#' angular-area coverage of a segmented field of view.
#'
#' @details
#' Every mask selects rows from the same `sampling_indices` table. Radiometric
#' variability is the median of `cv` at its selected indices, extracted without
#' window averaging. Coverage is the fraction of reference angular area in
#' segments containing at least one selected sample, as in
#' [assess_sampling_uniformity()]. The domain defined by `seg`, `w`, and optional
#' `m` is fixed across candidates. `bin_list` selects sky samples within that
#' domain; it does not change the reference areas. Areas are calculated once.
#'
#' Both criteria are independently min-max scaled across masks with at least
#' one sample. The Coverage-Consistency Metric is
#' \deqn{\mathrm{CCM}_i = (1 - \mathrm{coverage}_i)\,w +
#'   \mathrm{consistency}_i\,(1-w),}
#' where the two criteria are normalized and \eqn{w} is `coverage_weight`.
#'
#' The function selects the mask with the lowest CCM after ranking.
#' A constant criterion adds
#' the same contribution for every candidate and does not affect their ranking.
#' If all evaluable masks have the same score, the first is returned with a
#' warning. A single evaluable mask is returned without a warning.
#'
#' Masks selecting no samples have zero coverage, missing consistency and
#' missing CCM. They are reported with a warning and excluded from selection
#' and normalization. If none selects samples, the function stops. Nonfinite
#' CV values at selected indices cause an error.
#'
#' @note
#' Shrinking sky gaps can remove mixed pixels near their boundaries, reducing
#' local variability at the cost of losing samples and coverage. Low CV is the
#' chosen indicator of radiometric consistency; uniform non-sky regions can also
#' have low CV, but most binarization techniques classify them correctly so this
#' is not considered a design flaw.
#'
#' The candidate set determines the extrema used for normalization. Adding or
#' removing a candidate can change the selected mask even when the original
#' masks and their measurements remain unchanged. Inspect the original metrics
#' with `return = "all"` when comparing candidate sets.
#'
#' @param cv single-layer [terra::SpatRaster-class], typically from [cv_image()].
#'   Values at selected indices must be finite and nonnegative.
#' @param sampling_indices `data.frame` with positive integer-like `col` and
#'   `row` columns. Other columns are ignored. At least one row is required, but
#'   dense sampling is expected. All indices must lie within the reference
#'   domain defined by `seg`, `w`, and optional `m`. Repeated indices can be avoided via
#'   [rem_nearby_indices()], when provided are preserved  and contribute
#'   repeatedly to the median without, of course, increasing the number of
#'   covered cells.
#' @param bin_list nonempty list of single-layer logical, NA-free
#'   [terra::SpatRaster-class] objects. `TRUE` selects candidate sky pixels.
#'   List names identify candidates in the returned metrics.
#' @param coverage_weight numeric vector of length one in `[0, 1]`. Relative
#'   importance of coverage: zero selects by consistency alone, one by coverage
#'   alone.
#' @param return character vector of length one. `"mask"` returns the selected
#'   mask; `"all"` also returns its index and the metrics for every candidate.
#' @inheritParams assess_sampling_uniformity
#'
#' @return With `return = "mask"`, the selected element of `bin_list`, unchanged.
#'   With `return = "all"`, a list containing `mask`, `index` (its position in
#'   `bin_list`), and `data`, a data frame with:
#'   \describe{
#'     \item{`index`, `name`}{Candidate position and optional list name.}
#'     \item{`n_samples`}{Number of selected sampling rows.}
#'     \item{`consistency_raw`, `coverage_raw`}{Median CV and angular-area coverage
#'       before normalization.}
#'     \item{`consistency`, `coverage`}{Normalized criteria, or `NA` for masks
#'       without samples. Constant criteria are represented by zero.}
#'     \item{`ccm`}{Combined score; lower values are preferred.}
#'   }
#'
#' @section From gap boundaries to separated samples:
#' Use [distance_to_black()] to construct masks with progressively smaller sky
#' gaps. For an integer `k`, `distance_to_black(bin) > k` removes sky pixels
#' within Chebyshev distance `k` of black pixels; `k = 0` preserves the mask.
#' Candidates may originate from more than one initial binarization.
#'
#' Supply a common sampling table to evaluate the masks. Approximately uniform
#' angular sampling can be prepared with [fibonacci_directions()], followed by
#' [coord_from_internal_directions()] and [indices_from_coord()] for a radial geometry.
#' The sampling density and distribution affect both coverage and consistency.
#'
#' Once a mask is selected, calculate its distance to black and extract it at
#' the retained sampling indices. Pass that column to
#' [optimize_sampling_distance()] with `prioritize = "higher"`, then use the
#' same ranking in [rem_nearby_indices()]. This favors pixels deeper inside gaps
#' throughout the distance search and final selection. Distance to the boundary
#' is in pixels; separation between samples is angular. The resulting samples
#' can then be used to fit a sky radiance model.
#'
#' @seealso [distance_to_black()], [optimize_sampling_distance()],
#'   [rem_nearby_indices()], [cv_image()]
#'
#' @export
#'
#' @examples
#' # Synthetic image with two sky gaps and a smooth radiance gradient.
#' geometry <- list(model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' seg <- beckers_segmentation(d, 40)
#' w <- area_weight_image(g, geometry)
#' r <- terra::rast(g$off_axis)
#' xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
#' gaps <- (xy[, 1] - 35)^2 + (xy[, 2] - 35)^2 < 18^2 |
#'         (xy[, 1] - 52)^2 + (xy[, 2] - 17)^2 < 9^2
#' r[] <- ifelse(gaps, 100 + xy[, 1] / 5 + sin(xy[, 2] / 5), 10)
#' names(r) <- "Blue"
#' bin <- r > 50 & seg > 0
#' cv <- cv_image(r)
#'
#' distance <- distance_to_black(bin)
#' bin_list <- lapply(0:3, function(k) distance > k)
#' names(bin_list) <- paste0("shrink_", 0:3)
#' directions <- fibonacci_directions(3)
#' coord <- coord_from_internal_directions(directions, geometry)
#' sampling_indices <- indices_from_coord(coord$x, coord$y_down)
#' sampling_indices <- rem_unselected(sampling_indices, seg > 0)
#' sampling_indices <- rem_nearby_indices(sampling_indices, 1, space = "planar")
#'
#' choice <- optimize_sky_mask(cv, sampling_indices, seg, w, bin_list,
#'                             return = "all")
#' choice$data
#' plot(choice$mask)
#'
#' # Use the selected mask's boundary distance throughout angular thinning.
#' distance <- distance_to_black(choice$mask)
#' sampling_indices <- rem_unselected(sampling_indices, choice$mask)
#' sampling_indices <- add_directions(sampling_indices, d = d)
#' sampling_indices <- extract_dn(distance, sampling_indices, use_window = FALSE)
#' result <- optimize_sampling_distance(sampling_indices, seg, w,
#'   distances = seq(1, 45, 2), window_length = 3,
#'   col_id = "distance_to_black", prioritize = "higher")
#' if (!is.na(result$optimum)) {
#'   sky_indices <- rem_nearby_indices(sampling_indices, result$optimum,
#'     space = "spherical", col_id = "distance_to_black", prioritize = "higher")
#'   sky_samples <- extract_dn(r, sky_indices, use_window = FALSE)
#'   xy <- xy_from_indices(sky_samples$col, sky_samples$row, nrow(r))
#'   points(xy$x, xy$y, col = "violet", pch = 10)
#' }
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/DSCN4248.tif", package = "rcaiman")
#' metadata <- read.delim(gsub("tif", "tsv", path))
#' sensor_data <- read_caim(path)
#'
#' # Get specs
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' raw <- get_signature(my_entry, "raw")
#' cygm <- get_band_spec(raw, "cygm")
#' raw_geometry <- get_geometry_spec(cygm, "simple_method")
#' raw_radiometry <- get_radiometry_spec(raw_geometry, "photometric_sphere")
#'
#' # Represent image geometry
#' g <- set_internal_geometry(raw$dim, raw_geometry)
#' d <- direction_image(g)
#' m <- !is.na(g$off_axis)
#' w <- area_weight_image(g, raw_geometry)
#'
#' # Read and preprocess the image
#' caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
#' caim
#' r <- caim$Cyan
#' ff <- trend_surface_image(r, raw_radiometry$parameters[[metadata$aperture]]$Cyan, m = m)
#' r <- apply_flat_field(r, ff)
#'
#' # Binarize the image
#' query_directions <- fibonacci_directions(20)
#' thr <- apply_by_direction(r, g, query_directions, m, w = w,
#'                           fun = thr_isodata)
#' bin <- binarize_with_thr(r, thr$v1)
#'
#' # Prepare input arguments
#' seg <- beckers_segmentation(d, 200)
#' bin <- bin & seg > 0
#' cv <- cv_image(r)
#' distance <- distance_to_black(bin)
#' bin_list <- lapply(0:5, function(k) distance > k)
#' names(bin_list) <- paste0("shrink_", 0:5)
#'
#' directions <- fibonacci_directions(1/3)
#' directions <- select_directions_in_fov(directions, extract_fov(g))
#' sampling_indices <- coord_from_internal_directions(directions, raw_geometry) %>%
#'   add_indices()
#' sampling_indices <- rem_unselected(sampling_indices, seg > 0)
#' sampling_indices <- rem_nearby_indices(sampling_indices, 1, space = "planar")
#'
#' # Run the function and get the best sky mask
#' choice <- optimize_sky_mask(cv,
#'                             sampling_indices,
#'                             seg,
#'                             w,
#'                             bin_list,
#'                             coverage_weight = 0.5,
#'                             return = "all")
#' choice$data
#' bin <- choice$mask
#'
#' # Use the sky mask to sample
#' distance <- distance_to_black(bin)
#' sampling_indices <- rem_unselected(sampling_indices, bin)
#' sampling_indices <- add_directions(sampling_indices, d = d)
#' sampling_indices <- extract_dn(distance, sampling_indices, use_window = FALSE)
#'
#' result <- optimize_sampling_distance(
#'   sampling_indices,
#'   seg,
#'   w,
#'   distances = seq(1, 20, 0.5),
#'   window_length = 5,
#'   col_id = "distance_to_black",
#'   prioritize = "higher"
#' )
#'
#' if (!is.na(result$optimum)) {
#'   sky_indices <- rem_nearby_indices(
#'     sampling_indices,
#'     result$optimum,
#'     space = "spherical",
#'     col_id = "distance_to_black",
#'     prioritize = "higher"
#'   )
#'   sky_samples <- extract_dn(r, sky_indices, use_window = FALSE)
#' }
#' head(sky_samples)
#' xy <- sky_samples %>% add_coord() %>% add_xy(., nrow(r))
#'
#' # Visualize
#' plot(bin)
#' points(xy$x, xy$y, col = "violet", pch = 10)
#' }
optimize_sky_mask <- function(cv, sampling_indices, seg, w, bin_list,
                              coverage_weight = 0.5, return = "mask", m = NULL) {
  .assert_single_layer(cv)
  reference <- .sampling_area_reference(seg, w, m)
  .assert_same_geom(cv, seg)
  .check_vector(coverage_weight, "numeric", 1, sign = "nonnegative")
  if (coverage_weight > 1) stop("`coverage_weight` cannot be greater than one.", call. = FALSE)
  .assert_choice(return, c("mask", "all"))
  if (!is.list(bin_list) || length(bin_list) == 0) {
    stop("`bin_list` must be a nonempty list of masks.", call. = FALSE)
  }
  for (bin in bin_list) {
    .assert_logical_mask(bin)
    .assert_same_geom(bin, cv)
  }
  .check_columns(sampling_indices, required_columns = c("col", "row"))
  indices <- sampling_indices[, c("col", "row"), drop = FALSE]
  labels <- .sampling_reference_ids(indices, reference)
  if (nrow(indices) == 0) {
    stop("`sampling_indices` must contain at least one row.", call. = FALSE)
  }
  cv_values <- extract_dn(cv, indices, use_window = FALSE)[[3]]
  metrics <- lapply(bin_list, function(bin) {
    selected <- extract_dn(bin, indices, use_window = FALSE)[[3]]
    values <- cv_values[selected]
    .check_vector(values, "numeric", sign = "nonnegative", name = "selected CV values")
    list(n_samples = sum(selected),
         consistency_raw = if (length(values)) stats::median(values) else NA_real_,
         coverage_raw = sum(reference$q[unique(labels[selected])]))
  })
  candidate_names <- names(bin_list)
  if (is.null(candidate_names)) candidate_names <- rep(NA_character_, length(bin_list))
  data <- data.frame(index = seq_along(bin_list), name = candidate_names,
    n_samples = vapply(metrics, `[[`, numeric(1), "n_samples"),
    consistency_raw = vapply(metrics, `[[`, numeric(1), "consistency_raw"),
    coverage_raw = vapply(metrics, `[[`, numeric(1), "coverage_raw"))
  eligible <- data$n_samples > 0
  if (!any(eligible)) stop("No candidate mask selects sampling indices.", call. = FALSE)
  if (any(!eligible)) {
    warning("Masks without samples were excluded: ",
            paste(data$index[!eligible], collapse = ", "), ".", call. = FALSE)
  }
  normalize <- function(x) {
    out <- rep(NA_real_, length(x))
    values <- x[eligible]
    out[eligible] <- if (min(values) == max(values)) 0 else normalize_min_max(values)
    out
  }
  data$consistency <- normalize(data$consistency_raw)
  data$coverage <- normalize(data$coverage_raw)
  data$ccm <- (1 - data$coverage) * coverage_weight +
    data$consistency * (1 - coverage_weight)
  best <- which.min(data$ccm)
  if (sum(eligible) > 1 && length(unique(data$ccm[eligible])) == 1) {
    warning("All evaluable masks have the same score; selecting the first.", call. = FALSE)
  }
  if (return == "mask") bin_list[[best]] else
    list(mask = bin_list[[best]], index = best, data = data)
}
