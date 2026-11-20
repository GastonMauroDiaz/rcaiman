#' Interpolate values at viewing directions
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Interpolate sampled values using inverse angular distance weighting (IDW),
#' optionally blending them with a complementary source at the query directions.
#'
#' @details
#' [nn_spherical()] selects the neighbors and [apply_locally()] evaluates each
#' query. Both tables must use the same reference frame. Select external or
#' internal directions through `angle_columns`; a common rotation of all
#' directions preserves distances and interpolated values.
#'
#' A supplied `search_radius` defines a spherical cap centered on each query
#' (see [nn_spherical()]). The cap bounds the search region; `k` limits the
#' number of observations used from that region.
#'
#' Without active filling, weights are proportional to \eqn{1/\chi^p}, where
#' \eqn{\chi} is angular distance. If any retrieved neighbors coincide exactly
#' with a query, their mean value is returned and other neighbors are ignored.
#' With no eligible neighbors, the result is `NA`.
#'
#' @section Blending with a complementary source:
#' Supply a column of values at the query directions through `filling_source`.
#' These values must be in the same units as the sampled values. For example,
#' a model of relative radiance must be scaled to the sampled DN before blending.
#' The blend follows the form proposed by \insertCite{Lang2010;textual}{rcaiman}:
#' \deqn{B = \frac{\sum_j v_j b_j + W M}{\sum_j v_j + W},}
#' where \eqn{b_j} are neighboring observations, \eqn{M} is the complementary
#' value, and \eqn{W} is `filling_weight`.
#'
#' Here inverse-distance weights are rescaled using fixed angular bounds,
#' rather than the distances of the available neighbors:
#' \deqn{v_j = \operatorname{clip}_{[0,1]}
#'   \left(\frac{\chi_j^{-p} - \chi_{max}^{-p}}
#'   {\chi_{min}^{-p} - \chi_{max}^{-p}}\right).}
#' The bounds are `min_dist` and `search_radius`. A neighbor at or closer than
#' `min_dist` has weight one; a neighbor at `search_radius` has weight zero.
#' This preserves the influence of distance even when only one neighbor is
#' available. The implementation evaluates an equivalent expression that
#' avoids infinite weights at zero distance.
#'
#' A coincident observation therefore has finite weight in the blend and does
#' not override the complementary source. Without eligible neighbors, or when
#' all their rescaled weights are zero, the complementary value is returned.
#' With `rule = "all"`, this also applies when fewer than `k` neighbors are
#' available. A missing complementary value produces `NA` when filling is active.
#' A missing query direction always produces `NA`, even with a filling value.
#'
#' `filling_weight = 0` disables filling entirely and uses ordinary IDW, just
#' as `filling_source = NULL` does. `filling_weight` controls the complementary
#' source; it is not an observation weight or a pixel solid-angle weight.
#'
#' @param sampling_directions data frame containing angular columns selected
#'   by `angle_columns`, in radians, and a numeric column selected by `col_id`.
#'   Samples with missing values or directions are excluded before neighbor
#'   selection. Infinite values are not allowed. Additional columns are ignored.
#' @param query_directions data frame containing the same angular columns, in
#'   radians. Other columns are preserved. Must not already contain the column
#'   selected by `col_id` in `sampling_directions`.
#' @param col_id character or integer vector of length one. Name or position
#'   of the sampled-value column in `sampling_directions`.
#' @param filling_source optional character or integer vector of length one.
#'   Name or position of a numeric column in `query_directions` containing
#'   complementary values. Values must be finite or `NA`.
#' @param filling_weight nonnegative finite numeric vector of length one.
#'   Weight assigned to the complementary source. Zero disables filling.
#' @param p positive finite numeric vector of length one. IDW exponent.
#' @param min_dist positive finite numeric vector of length one. Lower angular
#'   bound in degrees for rescaling neighbor weights when filling is active.
#'   Must be smaller than `search_radius`. Ignored without active filling.
#' @param search_radius positive numeric vector of length one, at most 180,
#'   or `NULL`. Maximum angular distance in degrees for selecting neighbors.
#'   A finite radius is required when filling is active or `k = NULL`.
#' @inheritParams nn_spherical
#' @inheritParams apply_locally
#'
#' @return Data frame with the same rows, order, and columns as `query_directions`,
#'   plus the interpolated column, named as the sampled-value column. Queries
#'   without an estimate remain as `NA`; no rows are removed. Values are not
#'   normalized or clipped.
#'
#' @seealso [add_directions()], [nn_spherical()], [reproject_points()],
#'   [interpolate_planar()]
#' @export
#' @references \insertAllCited{}
#' @examples
#' sampling_directions <- sky_grid_centers(30)
#' sampling_directions$dn <- 100 + 20 * sampling_directions$zenith +
#'   5 * sin(sampling_directions$azimuth)
#' query_directions <- sky_grid_centers(15)
#' interpolated <- interpolate_spherical(
#'   sampling_directions, query_directions, col_id = "dn",
#'   k = 3, search_radius = 20
#' )
#' head(interpolated)
#'
#' # Complementary values for queries without nearby samples, or for blending.
#' query_directions$background <- 100 + 20 * query_directions$zenith
#' blended <- interpolate_spherical(
#'   sampling_directions, query_directions, col_id = "dn",
#'   filling_source = "background", filling_weight = 0.1,
#'   min_dist = 1, search_radius = 10
#' )
#' head(blended)
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
#' aperture <- as.character(metadata$aperture)
#' ff <- trend_surface_image(r, raw_radiometry$parameters[[aperture]]$Cyan, m = m)
#' r <- apply_flat_field(r, ff)
#'
#' # Read cie model
#' path <- system.file("external/cie_model.RDS", package = "rcaiman")
#' model <- readRDS(path)
#'
#' sky_cie <- cie_image(d, model$sun_angles, model$coef) * model$relative_radiance$zenith_dn
#' names(sky_cie) <- "above-canopy dn"
#' plot(sky_cie)
#'
#' # Binarize
#' bin <- binarize_with_thr(r/sky_cie, 0.5)
#' bin <- distance_to_black(bin) > 2
#' plot(bin)
#'
#' # Dense sampling of the sky DN
#' seg <- healpix_segmentation(d, 32)
#' seg[!bin] <- 0
#' features <- compute_features(r, seg, d, w)
#' sky_points <- features[features$effective_pixel_area > 3, ]
#' sky_points <- internal_directions_from_external(sky_points) %>%
#'   coord_from_internal_directions(., raw_geometry) %>% add_indices()
#' head(sky_points)
#' nrow(sky_points)
#'
#' # Prepare interpolation
#' fov <- extract_fov(g)
#' query_directions <- fibonacci_directions(1)
#' query_directions <- select_directions_in_fov(query_directions, fov)
#' query_directions <- external_directions_from_internal(query_directions)
#'
#' # Interpolate without blending
#' interpolated <- interpolate_spherical(
#'   sky_points,
#'   query_directions,
#'   col_id = "max",
#'   k = 3,
#'   search_radius = 20
#' )
#' head(interpolated)
#' interpolated <- coord_from_internal_directions(interpolated, raw_geometry)
#' plot(interpolate_planar(interpolated, "max", r,  search_radius = 10))
#'
#' # Get complementary values for blending
#' query_directions <- coord_from_internal_directions(query_directions, raw_geometry) %>%
#'   add_indices()
#' query_directions <- extract_dn(sky_cie, query_directions, use_window = FALSE)
#'
#' # Interpolate and blend
#' blended <- interpolate_spherical(
#'   sky_points,
#'   query_directions,
#'   col_id = "max",
#'   filling_source = "above-canopy dn",
#'   filling_weight = 0.01,
#'   min_dist = 1,
#'   rule = "all",
#'   search_radius = 20
#' )
#' head(blended)
#'
#' sky <- interpolate_planar(blended, "max", r,  search_radius = 10)
#' plot(sky)
#' plot(r/sky)
#' }
interpolate_spherical <- function(sampling_directions,
                                  query_directions,
                                  col_id,
                                  filling_source = NULL,
                                  filling_weight = 1,
                                  k = 3,
                                  p = 2,
                                  search_radius = 20,
                                  min_dist = 1,
                                  rule = "any",
                                  parallel = FALSE,
                                  cores = NULL,
                                  logical = TRUE,
                                  leave_free = 1,
                                  angle_columns = c("zenith", "azimuth")) {
  .check_columns(sampling_directions)
  if (is.numeric(col_id)) .check_vector(col_id, "integerish", 1, sign = "positive")
  col_id <- .resolve_required_column(col_id, sampling_directions)
  .check_columns(query_directions, reserved_columns = col_id)
  .check_vector(sampling_directions[[col_id]], "numeric", allow_na = TRUE,
                name = col_id)
  .check_vector(filling_weight, "numeric", 1, sign = "nonnegative")
  .check_vector(p, "numeric", 1, sign = "positive")
  if (!is.null(filling_source)) {
    if (is.numeric(filling_source)) {
      .check_vector(filling_source, "integerish", 1, sign = "positive")
    }
    filling_source <- .resolve_required_column(filling_source, query_directions)
    .check_vector(query_directions[[filling_source]], "numeric", allow_na = TRUE,
                  name = filling_source)
  }
  blend <- !is.null(filling_source) && filling_weight > 0
  if (blend) {
    .check_vector(min_dist, "numeric", 1, sign = "positive")
    .check_vector(search_radius, "numeric", 1, sign = "positive")
    if (min_dist >= search_radius) {
      stop("`min_dist` must be smaller than `search_radius` when filling is active.",
           call. = FALSE)
    }
  }
  # Validate the direction-column contract even for tables with no usable DN.
  .check_vector(angle_columns, "character", 2)
  if (anyDuplicated(angle_columns)) {
    stop("`angle_columns` must name two different columns.", call. = FALSE)
  }
  .spherical_neighbor_vectors(sampling_directions, angle_columns)
  valid_query <- stats::complete.cases(
    .spherical_neighbor_vectors(query_directions, angle_columns))
  samples <- sampling_directions[!is.na(sampling_directions[[col_id]]), , drop = FALSE]
  neighbors <- nn_spherical(samples, query_directions, k = k,
                             search_radius = search_radius, rule = rule,
                             parallel = parallel, cores = cores, logical = logical,
                             leave_free = leave_free, angle_columns = angle_columns)

  evaluate <- if (!blend) {
    function(data, distance, i) {
      if (!valid_query[i] || !nrow(data)) return(NA_real_)
      zero <- distance == 0
      if (any(zero)) return(mean(data[[col_id]][zero]))
      weights <- (min(distance) / distance)^p
      stats::weighted.mean(data[[col_id]], weights)
    }
  } else {
    lower <- .degree2radian(min_dist)
    upper <- .degree2radian(search_radius)
    function(data, distance, i) {
      if (!valid_query[i]) return(NA_real_)
      background <- query_directions[[filling_source]][i]
      if (is.na(background) || !nrow(data)) return(background)
      weights <- as.numeric(distance <= lower)
      between <- distance > lower & distance < upper
      chi <- distance[between]
      # Equivalent to fixed-bound min-max rescaling of inverse distances.
      weights[between] <- exp(p * log(lower / chi)) *
        (-expm1(p * log(chi / upper))) / (-expm1(p * log(lower / upper)))
      weights <- pmax(0, pmin(1, weights))
      if (any(!is.finite(weights))) {
        stop("Cannot compute finite blending weights from the supplied angular bounds.",
             call. = FALSE)
      }
      weights <- c(weights, filling_weight)
      stats::weighted.mean(c(data[[col_id]], background), weights / max(weights))
    }
  }
  # A zero-column query table retains rows without reserving names in the user's table.
  result <- apply_locally(samples, query_directions[, FALSE, drop = FALSE],
                           neighbors, fun = evaluate, parallel = parallel,
                           cores = cores, logical = logical, leave_free = leave_free)
  query_directions[[col_id]] <- result$output
  query_directions
}
