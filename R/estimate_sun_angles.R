#' Estimate the solar direction from sky brightness
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Estimate an initial solar direction from bright sky seen through canopy gaps.
#' Use localized brightness in veiled mode, or the spatial distribution of
#' bright segments in obscured mode.
#'
#' @param r single-layer numeric [terra::SpatRaster-class] of image brightness,
#'   after the appropriate radiometric preprocessing. Missing measurements are
#'   omitted. Measurements used as sky samples must otherwise be finite.
#' @param d external directions from [direction_image()], with `zenith` and
#'   `azimuth` in degrees, on the grid of `r`.
#' @param bin single-layer logical, NA-free [terra::SpatRaster-class] on the
#'   same grid. `TRUE` identifies sky samples; `FALSE` does not exclude a
#'   direction from the domain used to calculate cell areas and centers.
#' @param seg single-layer segmentation on the same grid, for example from
#'   [beckers_segmentation()] or [healpix_segmentation()]. Positive integer-like
#'   labels identify cells; zero and `NA` are excluded. IDs need not be
#'   consecutive. Use angularly compact cells; their resolution determines the
#'   spatial detail available to the estimator. Both modes use this argument.
#' @param w single-layer raster of nonnegative finite solid-angle weights,
#'   typically from [area_weight_image()]. Relative weights suffice. Zero and
#'   missing weights are excluded. All input rasters must share the same grid.
#' @param circumsolar_angle positive numeric vector of length one, at most
#'   180 deg. Search radius around the seed grouping in obscured mode. Groupings
#'   whose centers fall within this distance are retained. Ignored in veiled
#'   mode.
#' @param mode character vector of length one: `"obscured"` (default) or
#'   `"veiled"`.
#' @param m optional single-layer logical, NA-free raster restricting the domain
#'   to inspect. It is independent of `bin`. With `NULL`, use positive segment
#'   labels with available directions and positive weights above the horizon.
#' @param return character vector of length one. `"sun_angles"` returns the
#'   estimated direction; `"all"` also returns cell and grouping summaries.
#'
#' @details
#' Each cell's brightness is the 99th percentile of its available sky samples.
#' Vegetation pixels are omitted, not replaced by zero. The cell's area and
#' spherical centroid instead use its entire represented portion in the domain,
#' including directions hidden by vegetation. Cells without sky samples supply
#' no brightness evidence. A cell's centroid is not the centroid of its gaps.
#'
#' In **veiled mode**, the cell with the highest brightness supplies the solar
#' direction through its area-weighted spherical centroid. Ties select the
#' smallest cell ID. A fine segmentation permits a localized estimate when the
#' disk is visible, veiled, or covered by a sunscreen
#' \insertCite{Kuusk2007}{rcaiman}.
#'
#' In **obscured mode**, cells at or above the 90th percentile of cell brightness
#' are selected (one brightness value per sampled cell, including ties).
#' Adjacent bright cells form groupings using four-neighbor raster adjacency.
#' Each grouping is scored by `0.25 * area + 0.75 * brightness`, with both
#' criteria independently min-max scaled; a constant criterion contributes
#' equally to every grouping. Area is the sum of represented cell areas.
#' Grouping brightness is the area-weighted median of cell brightness values.
#' The highest-scoring grouping is the seed; ties select the grouping containing
#' the smallest cell ID. Retain groupings whose area-weighted spherical centers
#' are within `circumsolar_angle` of the seed center.
#'
#' The retained groupings are then separated into their constituent cells.
#' The estimate is the area-weighted spherical centroid of their represented
#' domain, including directions hidden by vegetation within those cells.
#' Uniform brightness provides no contrast for locating the sun: the result
#' then reflects the represented domain and the availability of sky samples.
#'
#' @note
#' This is an initial estimate for subsequent inspection and refinement with
#' [optimize_sun_angles()]. Bright cloud structures and incomplete sampling can
#' shift the selected region away from the sun. Its centroid describes the
#' observed bright region; it does not extrapolate a solar center below the
#' horizon. The example illustrates this limitation with a CIE sky whose solar
#' zenith angle is 92 deg.
#'
#' @return With `return = "sun_angles"`, a one-row data frame with `zenith` and
#'   `azimuth` in radians, as required by [cie_image()] and [fit_cie_model()].
#'   With `return = "all"`, a list with:
#'   \describe{
#'     \item{`sun_angles`, `method`}{Direction and method used: `brightest_cell`
#'       or `spherical_centroid`.}
#'     \item{`centroid`}{Area-weighted direction of the selected region.}
#'     \item{`cells`}{Table of cell IDs, relative areas, centroid directions in
#'       radians, sky-sample counts, brightness, grouping IDs, and flags
#'       `bright` and `selected`. Cells without sky samples have missing
#'       brightness.}
#'     \item{`groups`}{Obscured-mode grouping areas, brightness, centroid
#'       directions in radians, scores, and seed/selection flags; `NULL` in
#'       veiled mode.}
#'   }
#' @seealso [beckers_segmentation()], [healpix_segmentation()],
#'   [area_weight_image()], [optimize_sun_angles()]
#' @references \insertAllCited{}
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(160, 160),
#'   normalization_radius = 158,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(320, 320), geometry)
#' d <- direction_image(g)
#' w <- area_weight_image(g, geometry)
#' seg <- beckers_segmentation(d, 2400)
#' sun_angles <- data.frame(zenith = 92 * pi/180, azimuth = 0)
#' r <- cie_image(d, sun_angles, sky_coef = as.numeric(cie_table[10, 1:5]))
#' plot(r)
#'
#' # Synthetic gaps: retain randomly selected chessboard regions as sky.
#' # Angular cells retain their represented area behind vegetation.
#' cb <- chessboard(r, 10)
#' bin <- cb %in% sample(unique(cb[]), 300)
#' plot(bin)
#'
#' result <- estimate_sun_angles(r, d, bin, seg, w, return = "all")
#' result$sun_angles * 180 / pi
#' result$method
#' # This centroid is above the horizon, even though the modeled sun is below it.
#' # Changing the sampled gaps or segmentation can change the initial estimate.
estimate_sun_angles <- function(r,
                                d,
                                bin,
                                seg,
                                w,
                                circumsolar_angle = 30,
                                mode = "obscured",
                                m = NULL,
                                return = "sun_angles") {

  .assert_single_layer(r)
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  for (x in list(z, bin, seg, w)) .assert_same_geom(r, x)
  .assert_logical_mask(bin)
  .assert_choice(mode, c("obscured", "veiled"))
  .assert_choice(return, c("sun_angles", "all"))
  if (mode == "obscured") {
    .check_vector(circumsolar_angle, "numeric", 1, sign = "positive")
    if (circumsolar_angle > 180) stop("`circumsolar_angle` cannot exceed 180 deg.", call. = FALSE)
  }
  zenith <- terra::values(z, mat = FALSE)
  azimuth <- terra::values(a, mat = FALSE)
  # Azimuth is immaterial at the exact zenith.
  azimuth[!is.na(zenith) & zenith == 0 & is.na(azimuth)] <- 0
  domain <- !is.na(zenith) & zenith <= 90 & !is.na(azimuth)
  if (!is.null(m)) {
    .assert_logical_mask(m); .assert_same_geom(r, m)
    domain <- domain & terra::values(m, mat = FALSE)
  }
  reference <- .sampling_area_reference(seg, w,
    terra::setValues(terra::rast(r), domain))
  labels <- terra::values(reference$seg, mat = FALSE)
  domain <- !is.na(labels)
  weights <- terra::values(w, mat = FALSE)
  values <- terra::values(r, mat = FALSE)
  sky <- domain & terra::values(bin, mat = FALSE) & !is.na(values)
  .check_vector(values[sky], "numeric", name = "sky measurements")
  if (!any(sky)) stop("No sky measurements are available in the selected domain.", call. = FALSE)
  index <- match(labels, reference$ids)
  xyz <- .sun_unit_vectors(zenith[domain] * pi/180, azimuth[domain] * pi/180)
  area <- rowsum(weights[domain], index[domain], reorder = TRUE)[, 1]
  moments <- rowsum(xyz * weights[domain], index[domain], reorder = TRUE)
  norms <- sqrt(rowSums(moments^2))
  if (any(!is.finite(norms) | norms <= area * 1e-12)) {
    stop("A cell has an undefined spherical centroid.", call. = FALSE)
  }
  centers <- moments / norms
  directions <- .sun_vector_direction(centers)
  n <- length(area)
  sample_values <- split(values[sky], index[sky])
  brightness <- rep(NA_real_, n)
  brightness[as.integer(names(sample_values))] <- vapply(sample_values,
    stats::quantile, numeric(1), probs = .99, names = FALSE)
  cells <- data.frame(id = reference$ids, area = area, directions,
    n_samples = tabulate(index[sky], nbins = n), brightness = brightness,
    group = NA_integer_, bright = FALSE, selected = FALSE)
  groups <- NULL
  if (mode == "veiled") {
    best <- which.max(brightness)
    cells$bright[best] <- cells$selected[best] <- TRUE
    sun <- directions[best, , drop = FALSE]
    centroid <- sun
    method <- "brightest_cell"
  } else {
    threshold <- stats::quantile(brightness, .9, na.rm = TRUE, names = FALSE)
    cells$bright <- !is.na(brightness) & brightness >= threshold
    bright <- which(cells$bright)
    cell_group <- .sun_cell_groups(index, bright, nrow(r), ncol(r))
    cells$group[bright] <- cell_group
    group_area <- as.numeric(rowsum(area[bright], cell_group, reorder = TRUE))
    group_moments <- rowsum(moments[bright, , drop = FALSE], cell_group, reorder = TRUE)
    if (any(sqrt(rowSums(group_moments^2)) <= group_area * 1e-12)) {
      stop("A bright grouping has an undefined spherical centroid.", call. = FALSE)
    }
    group_direction <- .sun_vector_direction(group_moments)
    group_brightness <- vapply(split(bright, cell_group), function(i) {
      ord <- order(brightness[i]); i <- i[ord]
      brightness[i][which(cumsum(area[i] / sum(area[i])) >= .5)[1]]
    }, numeric(1))
    scale <- function(x) if (min(x) == max(x)) rep(1, length(x)) else normalize_min_max(x)
    score <- .25 * scale(group_area) + .75 * scale(group_brightness)
    seed <- which.max(score)
    distance <- calc_spherical_distance(group_direction$zenith, group_direction$azimuth,
      group_direction$zenith[seed], group_direction$azimuth[seed])
    groups <- data.frame(id = seq_along(group_area), area = group_area,
      brightness = unname(group_brightness), group_direction, score = unname(score),
      seed = seq_along(group_area) == seed,
      selected = distance <= circumsolar_angle * pi/180)
    cells$selected <- !is.na(cells$group) & cells$group %in% groups$id[groups$selected]
    selected <- which(cells$selected)
    total_moment <- colSums(moments[selected, , drop = FALSE])
    if (sqrt(sum(total_moment^2)) <= sum(area[selected]) * 1e-12) {
      stop("The selected region has an undefined spherical centroid.", call. = FALSE)
    }
    centroid <- .sun_vector_direction(matrix(total_moment, 1))
    sun <- centroid
    method <- "spherical_centroid"
  }
  rownames(sun) <- NULL; rownames(centroid) <- NULL
  if (return == "sun_angles") sun else list(sun_angles = sun, method = method,
    centroid = centroid, cells = cells, groups = groups)
}

.sun_unit_vectors <- function(zenith, azimuth) {
  cbind(sin(zenith) * cos(azimuth), sin(zenith) * sin(azimuth), cos(zenith))
}

.sun_vector_direction <- function(xyz) {
  xyz <- xyz / sqrt(rowSums(xyz^2))
  azimuth <- atan2(xyz[, 2], xyz[, 1]) %% (2*pi)
  azimuth[azimuth >= 2*pi] <- 0
  data.frame(zenith = atan2(sqrt(rowSums(xyz[, 1:2, drop = FALSE]^2)), xyz[, 3]),
             azimuth = azimuth)
}

# Connect cell IDs, so a masked split of one cell does not create two cells.
.sun_cell_groups <- function(index, bright, height, width) {
  x <- matrix(match(index, bright, nomatch = 0L), height, width, byrow = TRUE)
  edges <- rbind(cbind(as.vector(x[-height, , drop = FALSE]), as.vector(x[-1, , drop = FALSE])),
                 cbind(as.vector(x[, -width, drop = FALSE]), as.vector(x[, -1, drop = FALSE])))
  edges <- unique(edges[edges[, 1] > 0 & edges[, 2] > 0 & edges[, 1] != edges[, 2], , drop = FALSE])
  parent <- seq_along(bright)
  root <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_len(nrow(edges))) {
    a <- root(edges[i, 1]); b <- root(edges[i, 2])
    parent[max(a, b)] <- min(a, b)
  }
  roots <- vapply(seq_along(bright), root, integer(1))
  match(roots, sort(unique(roots)))
}
