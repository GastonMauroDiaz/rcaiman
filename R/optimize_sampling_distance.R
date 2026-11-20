#' Estimate optimal sampling distance
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Estimate an angular sampling distance by balancing uniformity and coverage
#' across a sequence of candidate distances.
#'
#' @details
#' For each value in `distances`, nearby sampling indices are removed using
#' [rem_nearby_indices()] in spherical space, and the retained subset is evaluated
#' with [assess_sampling_uniformity()]. Each distance starts from the same input
#' table and uses the same `col_id` and `prioritize` ranking. Without a ranking,
#' selection follows row order. Use the same table, ranking, and angular columns
#' when applying the estimated distance with [rem_nearby_indices()].
#'
#' When `window_length` is provided, both metric sequences are smoothed using
#' a centered moving average with equal weights. The window spans consecutive
#' evaluated distances. Boundary values are padded
#' by repeating edge values.
#'
#' The resulting Kullback--Leibler divergence (`kl`) and coverage sequences are
#' then independently scaled to \eqn{[0,1]} using [normalize_min_max()]. The
#' optimum is the first crossing of these normalized curves, estimated by linear
#' interpolation between adjacent evaluated distances.
#'
#' Apply the result with [rem_nearby_indices()] to inspect the
#' resulting sampling pattern.
#'
#' If either metric is constant its min-max scaling
#' is undefined. The function warns which metric did not vary and returns
#' `optimum = NA`. Original metrics are retained; the constant metric's
#' normalized column and `diff` contain `NA`.
#'
#' @note
#' The optimum expresses this curve-crossing criterion within the explored
#' interval. The interval both limits the search and determines the extrema used
#' to normalize the curves. Extending it can therefore move the crossing even
#' where the evaluated sampling patterns have not changed. Smoothing also
#' affects the crossing.
#'
#' As in [assess_sampling_uniformity()], `seg`, `w`, and optional `m` define
#' a fixed reference domain. Segment areas are calculated once and used for
#' every candidate distance, including regions that receive no samples.
#'
#' @param sampling_indices `data.frame` with `col` and `row` cell indices and
#'   the angular columns selected by `angle_columns`, in radians, as returned
#'   by [add_directions()]. All indices must select positive labels in
#'   `seg` within the reference domain defined by `w` and optional `m`.
#'   At least one row is required. Additional columns are
#'   ignored unless selected by `col_id`.
#' @param distances numeric vector. Candidate angular distances in degrees,
#'   strictly increasing and within (0, 90]. At least one value is required.
#' @param window_length optional numeric vector of length one. Positive odd
#'   integer defining the moving window length used to smooth the metrics.
#'
#' @inheritParams assess_sampling_uniformity
#' @inheritParams nn_spherical
#' @inheritParams rem_nearby_indices
#'
#' @section Prioritizing pixels within sky gaps:
#' After selecting a sky mask with [optimize_sky_mask()], extract
#' [distance_to_black()] at its sampling indices and use that column with
#' `prioritize = "higher"`. This favors candidates farther from gap boundaries
#' at every evaluated distance, as well as in the final selection. Distance to
#' black is measured in pixels; `distances` controls angular separation in
#' degrees. See [optimize_sky_mask()] for the complete workflow.
#'
#' @return List with the following components:
#' \describe{
#'   \item{`optimum`}{Estimated optimal distance (degrees), or `NA` when either
#'     metric is constant after optional smoothing.}
#'   \item{`data`}{`data.frame` with columns:
#'     \describe{
#'       \item{`distance`}{Evaluated distances (degrees).}
#'       \item{`kl_raw`, `coverage_raw`}{Original metrics, before smoothing
#'         and normalization.}
#'       \item{`kl`, `coverage`}{Metrics after optional smoothing and min-max
#'         normalization. A constant metric is represented by `NA`.}
#'       \item{`diff`}{Difference `kl - coverage`.}
#'     }}
#' }
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' seg <- beckers_segmentation(d, 20)
#' w <- area_weight_image(g, geometry)
#'
#' # An illustrative mask; it should be a selected sky mask.
#' bin <- seg > 0
#' distance <- distance_to_black(bin)
#' sampling_indices <- expand.grid(col = seq(1, 80, 3), row = seq(1, 60, 3))
#' sampling_indices <- rem_unselected(sampling_indices, bin)
#' sampling_indices <- add_directions(sampling_indices, d = d)
#' sampling_indices <- extract_dn(distance, sampling_indices, use_window = FALSE)
#' head(sampling_indices)
#'
#' res <- optimize_sampling_distance(
#'   sampling_indices,
#'   seg,
#'   w,
#'   distances = seq(1, 45, 2),
#'   window_length = 3,
#'   col_id = "distance_to_black",
#'   prioritize = "higher"
#' )
#'
#' with(res$data, {
#'   plot(distance, kl, type = "l", col = "blue",
#'        xlab = "Angular distance (deg)", ylab = "Normalized metric")
#'   lines(distance, coverage, col = "red")
#'   if (!is.na(res$optimum)) abline(v = res$optimum)
#' })
#' legend("right", legend = c("kl", "coverage"),
#'        lty = 1, col = c("blue", "red"))
#'
#' if (!is.na(res$optimum)) {
#'   retained_indices <- rem_nearby_indices(sampling_indices,
#'                                          min_dist = res$optimum, space = "spherical",
#'                                          col_id = "distance_to_black", prioritize = "higher")
#'   assess_sampling_uniformity(retained_indices, seg, w)
#' }
#'
#' if (!is.na(res$optimum)) {
#'   plot(bin)
#'   xy <- retained_indices %>% add_coord() %>% add_xy(., nrow(bin))
#'   points(xy$x, xy$y)
#' }
optimize_sampling_distance <- function(sampling_indices,
                                       seg,
                                       w,
                                       distances,
                                       window_length = NULL,
                                       angle_columns = c("zenith", "azimuth"),
                                       col_id = NULL,
                                       prioritize = NULL,
                                       m = NULL) {

  reference <- .sampling_area_reference(seg, w, m)
  .check_vector(distances, "numeric", sign = "positive")
  .check_vector(window_length, "odd_integerish", 1, sign = "positive", allow_null = TRUE)
  if (length(distances) == 0 || any(distances > 90)) {
    stop("`distances` must contain at least one value within (0, 90] degrees.",
         call. = FALSE)
  }
  if (is.unsorted(distances, strictly = TRUE)) {
    stop("`distances` must be strictly increasing.", call. = FALSE)
  }
  # Check all indices before thinning could conceal an invalid sample.
  .assess_sampling_reference(sampling_indices, reference)
  if (nrow(sampling_indices) == 0) {
    stop("`sampling_indices` must contain at least one row.", call. = FALSE)
  }

  metrics <- lapply(distances, function(distance) {
    sp <- rem_nearby_indices(sampling_indices, distance, space = "spherical",
                             angle_columns = angle_columns,
                             col_id = col_id, prioritize = prioritize)
    .assess_sampling_reference(sp, reference)
  })
  kl_raw <- vapply(metrics, function(x) x$kl, numeric(1))
  cov_raw <- vapply(metrics, function(x) x$coverage, numeric(1))

  .moving_average <- function(y) {
    if (is.null(window_length)) return(y)
    pad <- floor(window_length / 2)
    padded <- c(rep(y[1], pad), y, rep(y[length(y)], pad))
    averaged <- stats::filter(padded, rep(1 / window_length, window_length),
                              sides = 2)
    as.numeric(averaged[seq_along(y) + pad])
  }
  curves <- list(kl = .moving_average(kl_raw), coverage = .moving_average(cov_raw))
  constant <- vapply(curves, function(y) min(y) == max(y), logical(1))
  scaled <- lapply(curves, function(y) {
    if (min(y) == max(y)) rep(NA_real_, length(y)) else normalize_min_max(y)
  })
  f <- scaled$kl - scaled$coverage
  opt <- NA_real_
  if (any(constant)) {
    warning("Cannot estimate `optimum`: ", paste(names(curves)[constant], collapse = " and "),
            " did not vary after optional smoothing over the evaluated distances.",
            call. = FALSE)
  } else {
    opt <- distances[which.min(abs(f))]
    cross <- which(diff(sign(f)) != 0)
    if (length(cross) > 0) {
      i1 <- cross[1]
      i2 <- i1 + 1
      opt <- distances[i1] - f[i1] * (distances[i2] - distances[i1]) /
        (f[i2] - f[i1])
    }
  }

  list(optimum = opt, data = data.frame(
    distance = distances,
    kl_raw = kl_raw, coverage_raw = cov_raw,
    kl = scaled$kl, coverage = scaled$coverage, diff = f
  ))
}
