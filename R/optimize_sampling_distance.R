#' Estimate optimal sampling distance
#'
#' @description
#' Estimate the optimal angular distance by balancing sampling uniformity
#' metrics across a sequence of distances.
#'
#' @details
#' For each value in `distances`, nearby sky points are removed using
#' [rem_nearby_points()], and sampling uniformity is assessed with
#' [assess_sampling_uniformity()].
#'
#' The resulting Kullback–Leibler divergence (`kl`) and coverage sequences are
#' scaled to \eqn{[0,1]} using [normalize_min_max()].
#'
#' The optimal distance is defined as the value minimizing the absolute
#' difference between both metrics.
#'
#' If a sign change is detected in the difference, a linear interpolation
#' between adjacent points refines the estimate.
#'
#' When `window_length` is provided, both metric sequences are smoothed using
#' a centered moving average. Boundary values are padded by repeating edge
#' values, which may slightly distort the curves near the extremes of
#' `distances`.
#'
#' @param sky_points `data.frame` with columns `"row"`, `"col"`, `"zen"`, and
#'   `"azi"`, as in the output of [add_directions()].
#' @param distances numeric vector. Angular distances (deg), strictly increasing
#'   and within the range (0, 90].
#' @param window_length optional numeric vector of length one. Odd integer
#'   defining the moving window length used to smooth the metric sequences.
#'
#' @inheritParams assess_sampling_uniformity
#'
#' @return List with the following components:
#' \describe{
#'   \item{`optimum`}{Estimated optimal distance (deg).}
#'   \item{`data`}{`data.frame` with columns:
#'     \describe{
#'       \item{`distance`}{Evaluated distances (deg).}
#'       \item{`kl`}{Normalized Kullback–Leibler divergence.}
#'       \item{`coverage`}{Normalized coverage.}
#'       \item{`diff`}{Difference `kl - coverage`.}
#'     }}
#' }
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
#' sky_points <- sample_sky_points(r, bin, seg, dist_to_black = 3)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- add_directions(sky_points, z, a)
#' head(sky_points)
#'
#' distances <- seq(1, 20, by = 0.5)
#' res <- optimize_sampling_distance(sky_points, seg,
#'                                   distances,
#'                                   window_length = 5)
#'
#' # Plot results
#' with(res$data, {
#'   plot(distance, kl, type = "l", col = "blue",
#'        xlab = "Angular distance (deg)", ylab = "Normalized metric")
#'   lines(distance, coverage, col = "red")
#'   abline(v = res$optimum)
#' })
#'
#' legend("right", legend = c("kl", "coverage"),
#'        lty = 1, col = c("blue", "red"))
#' }
optimize_sampling_distance <- function(sky_points, equal_area_seg, distances,
                                       window_length = NULL) {

  .check_columns(sky_points, required_columns = c("row", "col", "zen", "azi"))
  .assert_equal_area_seg(equal_area_seg)
  .check_vector(distances, "numeric", sign = "positive")
  .check_vector(window_length, "odd_integerish", 1, sign = "positive", allow_null = TRUE)

  .moving_average <- function(y, window_length = 5) {
    if (!is.null(window_length)) {
      pad <- floor(window_length / 2)
      y_padded <- c(rep(y[1], pad), y, rep(y[length(y)], pad))
      y <- stats::filter(y_padded,
                         rep(1 / window_length, window_length), sides = 2)
      y <- y[!is.na(y)]
    }
    y
  }

  if (any(distances > 90)) {
    stop("`distances` must be within the range (0, 90] degrees.")
  }

  if (is.unsorted(distances, strictly = TRUE)) {
    stop("`distances` must be strictly increasing.")
  }

  metrics <- lapply(distances, function(d) {
    sp <- rem_nearby_points(sky_points, d, space = "spherical")
    assess_sampling_uniformity(sp, equal_area_seg)
  })

  kl <- vapply(metrics, function(x) x$kl, numeric(1))
  cov <- vapply(metrics, function(x) x$coverage, numeric(1))

  if (!is.null(window_length)) {
    kl <- .moving_average(kl, window_length)
    cov <- .moving_average(cov, window_length)
  }

  kl <- normalize_min_max(kl)
  cov <- normalize_min_max(cov)

  f <- kl - cov

  # discrete optimum
  idx <- which.min(abs(f))
  opt <- distances[idx]

  # refine with linear interpolation if possible
  cross <- which(diff(sign(f)) != 0)
  if (length(cross) > 0) {
    i1 <- cross[1]
    i2 <- i1 + 1

    x1 <- distances[i1]
    x2 <- distances[i2]

    y1 <- f[i1]
    y2 <- f[i2]

    opt <- x1 - y1 * (x2 - x1) / (y2 - y1)
  }

  list(
    optimum = opt,
    data = data.frame(
      distance = distances,
      kl = kl,
      coverage = cov,
      diff = f
    )
  )
}
