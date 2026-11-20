#' HEALPix segmentation
#'
#' @description
#' Assign observed external directions to cells of the HEALPix tessellation
#' on the sky hemisphere.
#'
#' @details
#' HEALPix (Hierarchical Equal Area isoLatitude Pixelization) partitions the
#' sphere into twelve equal-area base cells and recursively subdivides them
#' while preserving equal area and hierarchical relationships between cells
#' \insertCite{Gorski2005;textual}{rcaiman}.
#'
#' The implementation was developed from the original HEALPix description and
#' informed by the open-source Healpix.jl project
#' \insertCite{Tomasi2021;textual}{rcaiman}.
#'
#' This function assigns each pixel to a HEALPix cell according to its zenith
#' and azimuth coordinates and returns the corresponding cell identifiers using
#' the NESTED indexing scheme, numbered from 1. Both this function and
#' [adaptive_healpix_segmentation()] convert geographic azimuth to the
#' mathematical longitude used by the tessellation: `pi / 2 - azimuth`
#' in radians. This preserves the existing orientation of HEALPix cells.
#' Cell IDs do not depend on the image extent or acquisition orientation.
#'
#' Missing zenith or azimuth angles and directions below the horizon receive
#' zero. This includes the exact zenith when its azimuth is undefined; no
#' arbitrary polar cell is chosen.
#'
#' The horizon intersects some cells, so their upper-hemisphere portions have
#' less area than complete cells. The field of view and deliberate masks can
#' further restrict the represented portions, which need not have equal areas.
#'
#' @note
#' The proportion of horizon-intersected cells decreases as `nside`
#' increases. As with [beckers_segmentation()], use [area_weight_image()]
#' to account for the represented portions. [assess_sampling_uniformity()]
#' accepts both constructors and uses these areas as its sampling reference.
#'
#' @param nside integer vector of length one. Resolution parameter controlling
#'   the number of spherical regions. Larger values increase spatial
#'   subdivision. Powers of two (`1`, `2`, `4`, `8`, ...) follow the standard
#'   HEALPix NESTED indexing scheme. At most 8192, to keep cell IDs within
#'   the integer range used by the implementation.
#'
#' @inheritParams sky_grid_segmentation
#'
#' @return Integer [terra::SpatRaster-class] with one layer and the same
#'   geometry as `d`, with zero outside the segmentation. The attribute
#'   `cell_solid_angle` is `4 * pi / (12 * nside^2)`, in steradians, for a
#'   complete HEALPix cell. Observed portions of cells may cover less area.
#'
#' @seealso [beckers_segmentation()]
#'
#' @export
#'
#' @references
#' \insertAllCited{}
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
#' seg <- healpix_segmentation(d, nside = 4)
#' plot(seg)
healpix_segmentation <- function(d, nside = 1) {
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  .check_vector(nside, "integerish", 1, sign = "positive")
  if (!nside %in% 2^(0:13)) {
    stop("`nside` must be a power of two, at most 8192.", call. = FALSE)
  }
  zenith <- terra::values(z, mat = FALSE)
  azimuth <- terra::values(a, mat = FALSE)
  zenith[!is.na(zenith) & zenith > 90] <- NA_real_
  labels <- healpix_cpp(zenith, azimuth, nside)
  labels[is.na(labels)] <- 0L
  out <- terra::rast(z)
  terra::values(out) <- labels
  names(out) <- "HEALPix segmentation"
  attr(out, "cell_solid_angle") <- 4 * pi / (12 * nside^2)
  out
}
