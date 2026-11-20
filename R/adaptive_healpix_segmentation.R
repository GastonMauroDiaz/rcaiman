#' Adaptive HEALPix segmentation
#'
#' @description
#' Segment a hemispherical image into adaptive spherical regions by recursively
#' subdividing HEALPix cells according to local digital-number variability.
#'
#' @details
#' The procedure starts from an initial HEALPix tessellation
#' ([healpix_segmentation()]) and recursively subdivides regions according to
#' local radiometric variability. Regions with high variability tend to be
#' represented at finer spatial scales, whereas homogeneous regions remain
#' larger. The partition is defined in external directions from `d`, using
#' the same orientation as [healpix_segmentation()]. The split decisions and
#' returned labels depend on the observed data; labels are not fixed HEALPix IDs.
#'
#' Missing image values or directions, and directions below the horizon, are
#' excluded before fitting the tree and receive zero. An undefined azimuth
#' at the exact zenith is also excluded.
#'
#' HEALPix provides an equal-area spherical tessellation, reducing the
#' geometric distortions introduced by image-space partitioning when working
#' with hemispherical imagery. The recursive structure behaves similarly to a
#' quadtree, but subdivision occurs in spherical coordinates instead of image
#' coordinates.
#'
#' For each candidate split, a variability criterion is computed from the
#' standard deviation of pixel values within the parent and child regions:
#'
#' \deqn{
#' \Delta =
#' \sum_{i=1}^{4}\sigma_i -
#' \sigma_p
#' }
#'
#' where \eqn{\sigma_p} is the standard deviation of pixel values within the
#' parent region and \eqn{\sigma_i} is the standard deviation within child
#' region \eqn{i}.
#'
#' The `scale_parameter` determines the minimum increase in variability required
#' to justify additional subdivision. The seminal concept of this implementation
#' can be found in \insertCite{Diaz2015;textual}{rcaiman}.
#'
#' Recursive subdivision is accepted only when:
#'
#' \itemize{
#'   \item \eqn{\Delta > scale\_parameter}
#'   \item at least two child regions are non-empty
#'   \item child regions contain at least `min_size_px` pixels
#'   \item `max_splittings` has not been reached
#' }
#'
#' Increasing `scale_parameter` reduces the tendency toward additional
#' subdivision and generally produces larger regions, whereas smaller values
#' increase sensitivity to local variability and result in finer segmentations.
#'
#' The maximum subdivision depth is controlled by `max_splittings`, although
#' effective depth may be lower if regions become too small according to
#' `min_size_px`.
#'
#' @param r numeric [terra::SpatRaster-class] of one layer.
#' @param scale_parameter numeric vector of length one. Minimum increase in
#'   summed within-child standard deviation required to accept a subdivision.
#'   Larger values reduce splitting.
#' @param max_splittings numeric vector of length one. Must be positive and interpreted as
#'   integer. Maximum number of
#'   recursive subdivisions allowed, at most 13 to keep HEALPix IDs within
#'   the integer range used by the implementation.
#' @param min_size_px numeric vector of length one. Must be positive and interpreted as
#'   integer. Minimum number of pixels
#'   required in child regions for subdivision to be considered valid.
#'
#' @inheritParams sky_grid_segmentation
#'
#' @return Integer [terra::SpatRaster-class] with one layer and the same
#'   geometry as `r`, with zero outside the segmentation.
#'
#' @export
#'
#' @seealso [merge_segments()]
#'
#' @references \insertAllCited{}
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
#' r <- g$off_axis
#' seg <- adaptive_healpix_segmentation(r, d, scale_parameter = 0.1)
#' plot(seg)
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/DSCN4248.tif", package = "rcaiman")
#' sensor_data <- read_caim(path)
#'
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' raw <- get_signature(my_entry, "raw")
#' cygm <- get_band_spec(raw, "cygm")
#'
#' caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
#' caim
#' r <- caim$Cyan
#'
#' raw_geometry <- get_geometry_spec(cygm, "simple_method")
#' g <- set_internal_geometry(raw$dim, raw_geometry)
#' d <- direction_image(g)
#' seg <- adaptive_healpix_segmentation(r, d, scale_parameter = 25)
#' plot(seg)
#' plot(compute_mean(r, seg))
#' }
adaptive_healpix_segmentation <- function(r, d,
                        scale_parameter = 0,
                        max_splittings = 7,
                        min_size_px = 1) {
  .assert_single_layer(r)
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  .assert_same_geom(d, r)
  .check_vector(scale_parameter, "numeric", 1, sign = "nonnegative")
  .check_vector(max_splittings, "integerish", 1, sign = "positive")
  .check_vector(min_size_px, "integerish", 1, sign = "positive")
  if (max_splittings > 13) {
    stop("`max_splittings` must not exceed 13.", call. = FALSE)
  }
  zenith <- terra::values(z, mat = FALSE)
  azimuth <- terra::values(a, mat = FALSE)
  values <- terra::values(as.numeric(r), mat = FALSE)
  valid <- !is.na(zenith) & zenith <= 90 & !is.na(azimuth) & !is.na(values)
  .check_vector(values[valid], "numeric", name = "r at selected cells")
  zenith[!valid] <- NA_real_
  # The fixed HEALPix wrapper makes this same geographic-to-mathematical turn.
  longitude <- (pi / 2 - azimuth * pi / 180) %% (2 * pi)
  labels <- healpixtree_cpp(zenith * pi / 180, longitude, values,
    scale_parameter, max_splittings, min_size_px)
  labels[is.na(labels)] <- 0L
  seg <- terra::rast(r)
  terra::values(seg) <- labels
  names(seg) <- "healpix tree"
  seg
}
