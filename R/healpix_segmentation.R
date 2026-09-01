#' HEALPix segmentation
#'
#' @description
#' Partition a hemispherical image into equal-area spherical regions using the
#' HEALPix tessellation.
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
#' the standard NESTED indexing scheme.
#'
#' Cells belonging to the equatorial belt are intersected by the horizon and
#' contribute only half of their spherical area to the visible domain.
#' Consequently, regions close to the horizon have half the expected area within
#' the hemispherical image.
#'
#' @note
#' The proportion of horizon-intersected cells decreases as `nside`
#' increases.
#'
#' @param nside integer vector of length one. Resolution parameter controlling
#'   the number of spherical regions. Larger values increase spatial
#'   subdivision. Powers of two (`1`, `2`, `4`, `8`, ...) follow the standard
#'   HEALPix NESTED indexing scheme.
#'
#' @inheritParams add_directions
#'
#' @return Integer [terra::SpatRaster-class] with one layer and the same
#'   geometry as `z`. The object carries attribute `"cell_resolution_deg"`,
#'   calculated as `4 * pi / (12 * nside^2) * 180 / pi`.
#'
#' @seealso [equal_area_segmentation()]
#'
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#'
#' z <- zenith_image(600, equal_area_coef)
#' a <- azimuth_image(z)
#'
#' seg <- healpix_segmentation(
#'   z = z,
#'   a = a,
#'   nside = 1
#' )
#'
#' plot(seg)
#'
#' seg2 <- healpix_segmentation(
#'   z = z,
#'   a = a,
#'   nside = 2^6
#' )
#'
#' plot(seg2)
#'
#' }
#'
healpix_segmentation <- function(
    z,
    a,
    nside = 1
) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .is_power_of_two <- function(x) {
    x %in% unlist(lapply(0:100, function(x) 2^x))
  }
  .check_vector(nside, "integerish", 1, sign = "positive")
  if (nside > 2^30) {
    stop("`nside` not supported")
  }
  if (!.is_power_of_two(nside)) {
    stop("`nside` must be power of two")
  }

  z_values <- terra::values(z)
  a_values <- terra::values(a)

  labels <- healpix_cpp(
    z_values,
    a_values,
    nside
  )

  out <- z
  terra::values(out) <- labels
  out[is.na(out)] <- 0
  names(out) <- "HEALPix segmentation"
  attr(out, "cell_resolution_deg") <- 4*pi/(12*nside^2) * 180/pi

  out
}
