#' Work in a panoramic representation of the sky
#'
#' @description
#' Express external viewing directions as continuous panoramic coordinates,
#' process values on that grid, and transfer them to an image described by
#' external directions.
#'
#' @details
#' The unpadded domain covers zenith angles from 0 to 90 deg and geographic
#' azimuths from 0 to 360 deg. `x` increases from north toward east with
#' azimuth; `y_down` increases from zenith toward the horizon. One coordinate
#' unit spans `angle_width` degrees. The representation uses equal angular
#' steps, not equal solid angles.
#'
#' `to_panoramic()` transforms directions and creates a raster with one segment
#' per cell, also usable as a template for interpolation.
#' It does not aggregate or interpolate measurements. For example, use
#' [interpolate_planar()] to populate the template. Its search distances are
#' distances on this rectangular grid, not spherical distances. Interpolation
#' can extend beyond observed coverage; its extent is not a coverage mask.
#'
#' ### Aggregation by angular cells
#'
#' Use [aggregate_by_segment()] with the returned `raster` to summarize the
#' observations inside each panoramic cell. This
#' recovers the angular aggregation approach of the earlier `fisheye_to_pano()`
#' workflow, used in \insertCite{Diaz2021;textual}{rcaiman}. The summary may be
#' a mean, median, percentile, or another scalar function.
#'
#' To summarize the pixels of an image, supply every selected source pixel,
#' rather than a sparse sample. Each pixel then contributes once to its cell's
#' summary, without interpolation weights. Use `padding_deg = 0` because
#' aggregation does not need observations from neighboring cells. Cells without
#' observations remain `NA`. Assignment follows the containing-cell convention
#' of [aggregate_by_segment()], including its treatment of boundary coordinates.
#'
#' ### Padding
#'
#' The rectangular domain separates directions across the azimuth seam and
#' near zenith. Auxiliary copies of observations allow planar operations to
#' use data across these boundaries. Lateral copies shift `x` by a full turn.
#' To continue across zenith, a copy has negative zenith coordinate and its
#' azimuth is shifted by 180 deg. This combination represents the same
#' spherical direction. Lateral copies are added after the zenith copies,
#' including those needed at the upper corners. There is no horizon padding.
#'
#' Padding changes only the working coordinates. Source `zenith`, `azimuth`,
#' measured values, and other frame-independent columns remain attached to
#' each copy. These auxiliary rows are for panoramic processing, not new
#' observations. Padding does not remove the projection's metric distortion.
#'
#' ### Returning to an image
#'
#' Remove padding with `rem_panoramic_padding()` before calling
#' `from_panoramic()`. The latter selects the panoramic cell containing each
#' direction in `d`, without interpolation or angular aggregation. A direction
#' exactly at the horizon uses the last row. Directions below the horizon,
#' missing directions, and undefined azimuths at zenith produce `NA`; they
#' are not assigned values from an edge cell. Missing panoramic values remain
#' missing. This transfer does not undo previous resampling or processing.
#'
#' @param sampling_points data frame with numeric `zenith` and `azimuth`
#'   columns in radians, as returned by [add_directions()]. Zenith must be in
#'   `[0, pi]` and azimuth in `[0, 2 * pi)`. Missing angles are allowed.
#'   Original rows are retained; directions outside the sky hemisphere or
#'   with missing angles have `NA` working coordinates and generate no padding.
#'   Existing `x` and `y_down` are replaced, and frame-dependent columns are
#'   removed with [keep_coord()]. Other columns are preserved.
#' @param angle_width positive numeric vector of length one. Angular width
#'   of one panoramic cell in degrees. Must divide 90 into a whole number.
#' @param padding_deg integer-like numeric vector of length one in `[0, 90]`.
#'   Padding width in degrees. Must be a whole number of panoramic cells.
#'   Use zero for no padding.
#' @param pano numeric [terra::SpatRaster-class] with one or more layers,
#'   on the standardized grid from [set_rcaiman_geometry()]. Its dimensions
#'   must match `angle_width` and, for `rem_panoramic_padding()`, `padding_deg`.
#'   For `from_panoramic()`, supply an unpadded raster with cell values.
#' @param d [terra::SpatRaster-class] from [direction_image()], describing
#'   the destination, with `zenith` and `azimuth` layers in degrees, in that
#'   order. Its grid may be rectangular and its acquisition inclined.
#'
#' @return `to_panoramic()` returns a list containing:
#'   \describe{
#'     \item{`sampling_coord`}{A data frame with continuous `x`, `y_down`
#'       coordinates and preserved columns. Original rows appear first,
#'       followed by top, right-margin, and left-margin copies.}
#'     \item{`raster`}{A single-layer segmentation on the standardized image
#'       grid, with one segment per cell. Labels run from 1 to the number of
#'       cells, from left to right and top to bottom, including padding cells.
#'       It can also be used as a template for [interpolate_planar()], which
#'       ignores its cell values.}
#'   }
#'   `rem_panoramic_padding()` returns the central unpadded raster, with
#'   standardized geometry. `from_panoramic()` returns a raster on the grid
#'   of `d`. Both preserve the number and names of layers in `pano`.
#' @seealso [interpolate_planar()], [aggregate_by_segment()],
#'   [triangulate_planar()], [direction_image()]
#' @references \insertAllCited{}
#' @name panoramic
#' @rdname panoramic
#' @export
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
#' sampling_indices <- expand.grid(col = seq(1, ncol(d), 3),
#'                                 row = seq(1, nrow(d), 3))
#' m <- !is.na(d$zenith) & !is.na(d$azimuth) & d$zenith <= 90
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' samples <- add_directions(sampling_indices, d = d)
#' # Illustrative values; use extract_dn() for image measurements.
#' samples$dn <- cos(samples$zenith)
#' pano <- to_panoramic(samples, angle_width = 5, padding_deg = 10)
#' interpolated <- interpolate_planar(pano$sampling_coord, r = pano$raster,
#'   search_radius = 3, parallel = FALSE)
#' unpadded <- rem_panoramic_padding(interpolated, 5, 10)
#' result <- from_panoramic(unpadded, d, angle_width = 5)
#' plot(result)
#'
#' # Aggregate all selected image pixels instead of interpolating samples.
#' r <- cos(d$zenith * pi / 180)
#' names(r) <- "dn"
#' sampling_indices <- expand.grid(col = seq_len(ncol(d)), row = seq_len(nrow(d)))
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' samples <- add_directions(sampling_indices, d = d)
#' samples <- extract_dn(r, samples, use_window = FALSE)
#' pano <- to_panoramic(samples, angle_width = 10, padding_deg = 0)
#' aggregated <- aggregate_by_segment(pano$sampling_coord, pano$raster)
#' plot(aggregated)
#' # A percentile uses the same spatial assignment.
#' peaks <- aggregate_by_segment(pano$sampling_coord, pano$raster,
#'                               fun = quantile, probs = 0.99, na.rm = TRUE)
to_panoramic <- function(sampling_points, angle_width = 1, padding_deg = 15) {
  grid <- .panoramic_grid(angle_width, padding_deg)
  .check_columns(sampling_points, required_columns = c("zenith", "azimuth"))
  zenith <- sampling_points$zenith
  azimuth <- sampling_points$azimuth
  .check_vector(zenith, "numeric", allow_na = TRUE, name = "zenith")
  .check_vector(azimuth, "numeric", allow_na = TRUE, name = "azimuth")
  if (any(zenith < 0 | zenith > pi, na.rm = TRUE) ||
      any(azimuth < 0 | azimuth >= 2 * pi, na.rm = TRUE)) {
    stop("Directions must be in [0, pi] and [0, 2 * pi) radians.", call. = FALSE)
  }
  valid <- !is.na(zenith) & !is.na(azimuth) & zenith <= pi / 2
  coord <- sampling_points
  coord$x <- azimuth * 180 / pi
  coord$y_down <- zenith * 180 / pi
  coord$x[!valid] <- coord$y_down[!valid] <- NA_real_
  coord <- keep_coord(coord)

  top <- coord[which(valid & coord$y_down < padding_deg), , drop = FALSE]
  top$y_down <- -top$y_down
  top$x <- (top$x + 180) %% 360
  coord <- rbind(coord, top)
  right_margin <- coord[which(coord$x < padding_deg), , drop = FALSE]
  right_margin$x <- right_margin$x + 360
  left_margin <- coord[which(coord$x > 360 - padding_deg), , drop = FALSE]
  left_margin$x <- left_margin$x - 360
  coord <- rbind(coord, right_margin, left_margin)
  coord$x <- (coord$x + padding_deg) / angle_width
  coord$y_down <- (coord$y_down + padding_deg) / angle_width

  r <- terra::rast(nrows = grid$height, ncols = grid$width)
  r <- set_rcaiman_geometry(r)
  terra::values(r) <- seq_len(terra::ncell(r))
  list(sampling_coord = coord, raster = r)
}

#' @rdname panoramic
#' @export
rem_panoramic_padding <- function(pano, angle_width = 1, padding_deg = 15) {
  grid <- .panoramic_grid(angle_width, padding_deg)
  .check_panoramic_raster(pano, grid)
  e <- terra::ext(grid$margin, grid$margin + grid$columns,
                  0, grid$rows)
  set_rcaiman_geometry(terra::crop(pano, e))
}

#' @rdname panoramic
#' @export
from_panoramic <- function(pano, d, angle_width = 1) {
  grid <- .panoramic_grid(angle_width, 0)
  .check_panoramic_raster(pano, grid)
  if (!terra::hasValues(pano)) {
    stop("`pano` must contain cell values.", call. = FALSE)
  }
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  zenith <- terra::values(z, mat = FALSE)
  azimuth <- terra::values(a, mat = FALSE)
  valid <- which(!is.na(zenith) & !is.na(azimuth) & zenith <= 90)
  indices <- indices_from_coord(azimuth[valid] / angle_width,
                                zenith[valid] / angle_width)
  # Include the closed horizon boundary; all other out-of-domain rays stay NA.
  indices$row <- pmin(indices$row, grid$rows)
  # Guard rounding to the upper edge for azimuths immediately below 360 deg.
  indices$col <- pmin(indices$col, grid$columns)
  cells <- terra::cellFromRowCol(pano, indices$row, indices$col)
  values <- matrix(NA_real_, nrow = terra::ncell(d), ncol = terra::nlyr(pano))
  values[valid, ] <- terra::values(pano, mat = TRUE)[cells, , drop = FALSE]
  out <- terra::rast(z, nlyrs = terra::nlyr(pano))
  terra::values(out) <- values
  names(out) <- names(pano)
  out
}

.panoramic_grid <- function(angle_width, padding_deg) {
  rows <- .sky_interval_count(angle_width, 90)
  .check_vector(padding_deg, "integerish", 1, sign = "nonnegative")
  if (padding_deg > 90) stop("`padding_deg` must not exceed 90.", call. = FALSE)
  margin <- padding_deg / angle_width
  if (abs(margin - round(margin)) > 64 * .Machine$double.eps * max(1, margin)) {
    stop("`padding_deg` must span a whole number of panoramic cells.", call. = FALSE)
  }
  margin <- round(margin)
  list(rows = rows, columns = 4 * rows, margin = margin,
       height = rows + margin, width = 4 * rows + 2 * margin)
}

.check_panoramic_raster <- function(pano, grid) {
  if (!inherits(pano, "SpatRaster")) {
    stop("`pano` must be a SpatRaster.", call. = FALSE)
  }
  expected <- c(0, grid$width, 0, grid$height)
  if (nrow(pano) != grid$height || ncol(pano) != grid$width ||
      !isTRUE(all.equal(as.vector(terra::ext(pano)), expected,
                        check.attributes = FALSE))) {
    stop("`pano` must have the standardized extent and dimensions implied by ",
         "`angle_width` and padding.", call. = FALSE)
  }
  invisible(TRUE)
}
