#' Apply methods within angular windows
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Evaluate a function within overlapping angular windows and combine their
#' estimates on the original image grid.
#'
#' @param r numeric [terra::SpatRaster-class] with one layer. Missing image
#'   values do not supply samples.
#' @param g [terra::SpatRaster-class] from
#'   [set_internal_geometry()], with layers `off_axis` and `internal_azimuth`
#'   in degrees, on the same grid as `r`.
#' @param query_directions data frame with `off_axis` and `internal_azimuth`
#'   in radians, following [add_directions()]. Each row is the center of an
#'   angular window. For example, use [fibonacci_directions()] to generate them.
#' @param m logical, NA-free [terra::SpatRaster-class] on the same grid as `r`.
#'   Only pixels selected by this mask supply samples and receive estimates.
#' @param w optional single-layer raster of nonnegative, finite weights, or
#'   `NA`, on the same grid as `r`. For solid-angle weights, use
#'   [area_weight_image()]. Missing weights exclude the corresponding pixels;
#'   zero weights are passed to `fun`. With `NULL`, every supplied pixel has
#'   weight one.
#' @param fov positive numeric vector of length one, at most 360. Angular
#'   diameter of each window, in degrees. Its search radius is `fov / 2`.
#' @param laxity positive numeric vector of length one. Maximum deviation from
#'   the median, in scaled median absolute deviations, retained when combining
#'   estimates. Defaults to 2.5.
#' @param method character vector of length one. Either `"vector"` or
#'   `"single_layer_image"`.
#' @param fun function computing one or more estimates. With `method = "vector"`,
#'   it receives values and their corresponding weights as `function(x, w)`.
#'   With `method = "single_layer_image"`, it receives a local image, a logical
#'   mask, and a weight raster on the same local grid, as `function(r, m, w)`.
#'   Return a numeric scalar or vector, optionally named. Use `NA` for an
#'   unavailable estimate, or `NULL` to skip the window.
#'
#' @details
#' Each angular window is a spherical cap (see [nn_spherical()]) centered on
#' its query direction, with angular radius `fov / 2`. Its projection onto the
#' image need not be circular. The image coverage and `m` determine which
#' pixels are available within the cap.
#'
#' All image pixels with values, weights, and internal directions selected by `m` are
#' considered. [nn_spherical()] finds those within `fov / 2` degrees of each
#' query, with no limit on the number of neighbors. Missing directions are
#' excluded. Query directions remain continuous: they are not snapped to image
#' cells. A center outside the observed field of view can still define a window
#' containing selected pixels. Empty windows are skipped without calling `fun`.
#'
#' Selection uses angular distances in the internal frame; acquisition
#' orientation is not needed. Window centers and image directions must share
#' that frame. Partial windows use the available pixels without extrapolation
#' or a minimum coverage requirement. Weights are passed unchanged: `fun`
#' determines whether and how to use them. The function does not multiply
#' image values by weights or normalize weights within a window.
#'
#' In vector mode, values and weights are supplied in matching image cell
#' order. In image mode,
#' the local raster is the smallest rectangle containing the selected pixels,
#' with the image-coordinate convention of [set_rcaiman_geometry()]. Values
#' outside the window or mask are `NA`; the local mask selects the supplied
#' pixels. The weight raster has the same `NA` pattern as the local image.
#' The original pixel arrangement is preserved without reprojection.
#'
#' Each estimate is assigned to every sampled pixel in its window. Where
#' windows overlap, estimates are combined separately for each output variable.
#' Their median and median absolute deviation (MAD, scaled by 1.4826) define
#' the retained values, which are averaged. If MAD is zero, all values are
#' averaged. Weights supplied through `w` do not change this combination.
#' Counts record nonmissing contributions before filtering, not sums of weights.
#'
#' Returned vectors must have consistent lengths and names across windows;
#' unnamed elements receive names `v1`, `v2`, and so on. An unnamed scalar `NA`,
#' `NULL`, or an empty vector skips the whole window. A named or multivalue
#' result can contain `NA` for individual estimators. Errors in `fun` stop the
#' calculation with the query row number and the original error message.
#'
#' @section Methods that can fail in individual windows:
#' A method may not produce an estimate for every angular window. To skip
#' windows where the chosen method raises an error, handle that error inside
#' `fun`. For example, extract the middle threshold from [thr_two_corner()]
#' and return `NA` when it fails.:
#'
#' ```r
#' .fun <- function(x, w) {
#'   tryCatch({
#'     h <- build_histogram(x, w, n_bins = 256)
#'     peaks <- find_macfarlane_peaks(h)
#'     thr_two_corner(smooth_histogram(h), peaks)$tm
#'   }, error = function(e) NA)
#' }
#' thr <- apply_by_direction(
#'   r,
#'   g,
#'   fibonacci_directions(20),
#'   !is.na(g$off_axis),
#'   fun = .fun,
#'   w = w,
#'   fov = 40
#' )
#' ```
#'
#' Other overlapping windows can still supply estimates to those pixels.
#' The layer `n` records the nonmissing contributions before MAD filtering.
#'
#' The *Examples* section shows the use of fallback methods.
#'
#' @return A numeric [terra::SpatRaster-class] on the grid of `r`. Estimate
#'   layers retain the order of the first result, followed by count layers.
#'   A single estimator has a count layer named `n`; multiple estimators each
#'   have a count layer named `n_<estimator>`. Estimator names must not collide
#'   with count names. Pixels without estimates have `NA` values and zero
#'   counts, except that counts remain positive if MAD filtering rejects all
#'   their contributions. If no window returns a result defining the output
#'   variables, the output contains `v1` filled with `NA` and `n` with zeros.
#'
#' @note
#' This function is part of a manuscript currently under preparation.
#'
#' @seealso [fibonacci_directions()], [nn_spherical()], [binarize_with_thr()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(30, 20),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(60, 40), geometry)
#' w <- area_weight_image(g, geometry)
#' r <- g$off_axis + 20
#' names(r) <- "dn"
#' m <- !is.na(g$off_axis)
#' query_directions <- fibonacci_directions(40)
#'
#' estimates <- apply_by_direction(
#'   r, g, query_directions, m,
#'   fun = function(x, w) c(mean = mean(x), area_mean = weighted.mean(x, w)),
#'   w = w,
#'   method = "vector"
#' )
#' estimates
#' local_mean <- function(r, m, w) {
#'   selected <- terra::values(m, mat = FALSE)
#'   weighted.mean(terra::values(r, mat = FALSE)[selected],
#'                 terra::values(w, mat = FALSE)[selected])
#' }
#'
#' estimates <- apply_by_direction(
#'   r, g, query_directions, m,
#'   fun = local_mean,
#'   w = w,
#'   method = "single_layer_image"
#' )
#' estimates
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
#' w <- area_weight_image(g, raw_geometry)
#'
#' # Read and preprocess the image
#' caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
#' caim
#' r <- caim$Cyan
#' aperture <- as.character(metadata$aperture)
#' ff <- trend_surface_image(r,
#'                           raw_radiometry$parameters[[aperture]]$Cyan,
#'                           m = !is.na(g$off_axis))
#' r <- apply_flat_field(r, ff)
#'
#'
#' .fun <- function(r, m, ...) {
#'   r <- normalize_min_max(r) * 255
#'   r[!m] <- 0
#'   tryCatch(thr_nbhu(r,
#'                     support_weight = 1,
#'                     smooth_curve = TRUE, smooth_width = 7),
#'            error = function(e) NA)
#'
#' }
#' thr <- apply_by_direction(
#'   r,
#'   g,
#'   fibonacci_directions(20),
#'   !is.na(g$off_axis),
#'   fun = .fun,
#'   w = NULL,
#'   fov = 60,
#'   method = "single_layer_image"
#' )
#' plot(binarize_with_thr(r, thr$v1))
#'
#' # Weights should be used with radiometrically corrected images, otherwise
#' # the link between solid angle and response will produce uncertaining results.
#'
#' thr <- apply_by_direction(
#'   r,
#'   g,
#'   fibonacci_directions(20),
#'   !is.na(g$off_axis),
#'   fun = thr_isodata,
#'   w = w,
#'   fov = 60
#' )
#' plot(binarize_with_thr(r, thr$v1))
#'
#' .fun <- function(x, w) {
#'   tryCatch({
#'     h <- build_histogram(x, w, n_bins = 256)
#'     sigma <- 0
#'     repeat {
#'       sigma <- sigma + 1
#'       h <- smooth_histogram(h, sigma)
#'       peaks <- find_histogram_peaks(h)
#'       if (length(peaks) <= 2)
#'         break
#'     }
#'     thr_two_corner(h, peaks, diagnose = FALSE)$ts
#'   }, error = function(e) {
#'     if (length(peaks) != 2)
#'       return(NA)
#'     if (which.max(h$density[peaks]) == 1) {
#'       # thr <- h$mids[which.min(h$density[peaks[1]:peaks[2]])]
#'       tryCatch({
#'         h$density[1:which.min(h$density[peaks[1]:peaks[2]])] <- 0
#'         thr_rosin(h, peaks[1], tail = "right", diagnose = FALSE)
#'       }, error = function(e) NA
#'       )
#'     } else {
#'       # thr <- h$mids[which.min(h$density[peaks[1]:peaks[2]])]
#'       tryCatch({
#'         thr_rosin(h, peaks[2], tail = "right", diagnose = FALSE)
#'       }, error = function(e)
#'         NA)
#'     }
#'   })
#' }
#'
#' thr <- apply_by_direction(
#'   r,
#'   g,
#'   fibonacci_directions(5),
#'   !is.na(g$off_axis),
#'   fun = .fun,
#'   w = w,
#'   fov = 20
#' )
#' plot(binarize_with_thr(r, thr$v1))
#' }
apply_by_direction <- function(r, g, query_directions, m,
                               fun, w = NULL, fov = 40, laxity = 2.5,
                               method = "vector") {
  .assert_single_layer(r)
  .assert_spatraster(g)
  if (!identical(names(g), c("off_axis", "internal_azimuth"))) {
    stop("`g` must have layers `off_axis` and ",
         "`internal_azimuth`, in that order.", call. = FALSE)
  }
  .assert_same_geom(r, g)
  .assert_logical_mask(m)
  .assert_same_geom(r, m)
  if (!is.null(w)) {
    .assert_area_weight_image(w)
    .assert_same_geom(r, w)
  }
  .check_vector(fov, "numeric", 1, sign = "positive")
  if (fov > 360) stop("`fov` must not exceed 360 degrees.", call. = FALSE)
  .check_vector(laxity, "numeric", 1, sign = "positive")
  .assert_choice(method, c("vector", "single_layer_image"))
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)

  angle_columns <- c("off_axis", "internal_azimuth")
  # Check the full angular raster, including regions excluded by the mask.
  directions <- as.data.frame(terra::values(g)) * pi / 180
  xyz <- .spherical_neighbor_vectors(directions, angle_columns)
  values <- terra::values(r, mat = FALSE)
  weights <- if (is.null(w)) rep(1, terra::ncell(r)) else terra::values(w, mat = FALSE)
  selected <- terra::values(m, mat = FALSE) & !is.na(values) & !is.na(weights) &
    stats::complete.cases(xyz)
  cells <- which(selected)
  .check_vector(as.numeric(values[cells]), "numeric", name = "selected image values")
  neighbors <- nn_spherical(
    directions[cells, , drop = FALSE], query_directions,
    k = NULL, search_radius = fov / 2, angle_columns = angle_columns,
    parallel = FALSE
  )
  # Keep only row positions; distances are not needed for estimation.
  neighbors <- lapply(neighbors, function(x) sort(x$index))
  results <- vector("list", length(neighbors))
  estimator_names <- NULL
  for (i in seq_along(neighbors)) {
    hit <- cells[neighbors[[i]]]
    if (length(hit) == 0L) next
    result <- tryCatch({
      if (method == "vector") {
        fun(as.numeric(values[hit]), as.numeric(weights[hit]))
      } else {
        rc <- terra::rowColFromCell(r, hit)
        local <- terra::rast(nrows = diff(range(rc[, 1])) + 1,
                             ncols = diff(range(rc[, 2])) + 1)
        local <- set_rcaiman_geometry(local)
        local_values <- rep(NA_real_, terra::ncell(local))
        local_cells <- terra::cellFromRowCol(
          local, rc[, 1] - min(rc[, 1]) + 1, rc[, 2] - min(rc[, 2]) + 1
        )
        local_values[local_cells] <- values[hit]
        local <- terra::setValues(local, local_values)
        names(local) <- names(r)
        local_weights <- rep(NA_real_, terra::ncell(local))
        local_weights[local_cells] <- weights[hit]
        local_weights <- terra::setValues(terra::rast(local), local_weights)
        names(local_weights) <- "w"
        fun(local, !is.na(local), local_weights)
      }
    }, error = function(e) {
      stop("`fun` failed at query_directions row ", i, ": ",
           conditionMessage(e), call. = FALSE)
    })
    if (is.null(result) || length(result) == 0L ||
        (is.atomic(result) && length(result) == 1L &&
         is.null(names(result)) && is.na(result))) next
    if (is.logical(result) && all(is.na(result))) storage.mode(result) <- "double"
    if (!is.numeric(result) || is.complex(result) || !is.null(dim(result)) ||
        any(is.infinite(result))) {
      stop("`fun` must return a numeric vector with finite or NA values ",
           "(query_directions row ", i, ").", call. = FALSE)
    }
    if (is.null(names(result))) names(result) <- paste0("v", seq_along(result))
    if (anyNA(names(result)) || any(!nzchar(names(result))) ||
        anyDuplicated(names(result))) {
      stop("`fun` must return unique, nonempty estimator names ",
           "(query_directions row ", i, ").", call. = FALSE)
    }
    if (is.null(estimator_names)) {
      estimator_names <- names(result)
      count_names <- if (length(result) == 1L) "n" else paste0("n_", estimator_names)
      if (anyDuplicated(c(estimator_names, count_names))) {
        stop("Estimator names must not collide with count layer names.", call. = FALSE)
      }
    } else if (!setequal(names(result), estimator_names)) {
      stop("`fun` returned inconsistent estimators at query_directions row ",
           i, ".", call. = FALSE)
    }
    results[[i]] <- result[estimator_names]
  }
  if (is.null(estimator_names)) estimator_names <- "v1"
  count_names <- if (length(estimator_names) == 1L) "n" else paste0("n_", estimator_names)
  estimates <- counts <- vector("list", length(estimator_names))
  for (j in seq_along(estimator_names)) {
    valid <- which(vapply(results, function(x) length(x) > 0L && !is.na(x[j]), logical(1)))
    estimate_values <- rep(NA_real_, terra::ncell(r))
    count_values <- integer(terra::ncell(r))
    if (length(valid)) {
      lengths <- vapply(neighbors[valid], length, integer(1))
      contributing_cells <- cells[unlist(neighbors[valid], use.names = FALSE)]
      contributing_values <- rep(vapply(results[valid], function(x) x[j], numeric(1)),
                                 times = lengths)
      combined <- robust_grouped_mean_cpp(contributing_cells, contributing_values,
                                          laxity = laxity)
      estimate_values[combined$cell] <- combined$value
      count_values <- tabulate(contributing_cells, nbins = terra::ncell(r))
    }
    estimates[[j]] <- terra::setValues(terra::rast(r), estimate_values)
    counts[[j]] <- terra::setValues(terra::rast(r), count_values)
  }
  output <- terra::rast(c(estimates, counts))
  names(output) <- c(estimator_names, count_names)
  output
}
