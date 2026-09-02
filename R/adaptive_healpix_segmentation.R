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
#' larger.
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
#'   recursive subdivisions allowed.
#' @param min_size_px numeric vector of length one. Must be positive and interpreted as
#'   integer. Minimum number of pixels
#'   required in child regions for subdivision to be considered valid.
#'
#' @inheritParams add_directions
#'
#' @return Integer [terra::SpatRaster-class] with one layer and the same
#'   geometry as `r`.
#'
#' @export
#'
#' @seealso [merge_segments()]
#'
#' @references \insertAllCited{}
#'
#' @examples
#'
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
#' seg_tree <- adaptive_healpix_segmentation(
#'   r,
#'   z,
#'   a,
#'   scale_parameter = 5,
#'   max_splittings = 10,
#'   min_size_px = 5
#' )
#' plot(seg_tree)
#' }
adaptive_healpix_segmentation <- function(r, z, a,
                        scale_parameter = 0,
                        max_splittings = 7,
                        min_size_px = 1) {

  .assert_single_layer(r, allow_na = FALSE)
  .assert_zenith_image(z)
  .assert_same_geom(z, r)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(scale_parameter, "numeric", 1, sign = "nonnegative")
  .check_vector(max_splittings, "integerish", 1, sign = "positive")
  .check_vector(min_size_px, "integerish", 1, sign = "positive")

  # Las funciones C++ de HEALPix esperan radianes
  z_rad <- terra::values(z) * pi / 180
  a_rad <- terra::values(a) * pi / 180
  r_val <- terra::values(r)

  labels <- healpixtree_cpp(
    z = z_rad,
    a = a_rad,
    values = r_val,
    scale_parameter = scale_parameter,
    max_splittings = max_splittings,
    min_size_px = min_size_px
  )

  seg <- r
  terra::values(seg) <- labels
  names(seg) <- "healpix tree"
  seg[is.na(seg)] <- 0

  return(seg)
}
