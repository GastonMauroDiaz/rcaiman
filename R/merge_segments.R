.compute_adjacency <- function(seg,
                              w,
                              directed = FALSE,
                              ignore_label_0 = TRUE,
                              return_border_length = TRUE) {

  x <- terra::as.matrix(seg, wide = TRUE)
  w <- terra::as.matrix(w, wide = TRUE)

  out <- compute_adjacency_cpp(
    x = x,
    w = w,
    directed = directed,
    ignore_label_0 = ignore_label_0,
    return_border_length = return_border_length
  )

  rownames(out) <- NULL

  out
}

#' Multiresolution region merging
#'
#' @description
#' Iteratively merge adjacent image regions according to a weighted radiometric
#' heterogeneity criterion inspired by the multiresolution segmentation strategy
#' popularized by eCognition \insertCite{baatz2000}{rcaiman}.
#'
#' @details
#' This function implements a bottom-up region-merging procedure operating on
#' an existing labeled segmentation raster.
#'
#' The algorithm starts from an initial oversegmentation and iteratively merges
#' neighboring regions that minimize the heterogeneity of the resulting merged
#' object.
#'
#' Unlike classical raster-based implementations where each pixel contributes
#' equally, this implementation supports weighted radiometric statistics via
#' an area-weight image (`w`). This enables segmentation directly in any known
#' hemispherical projection.
#'
#' @section Multiresolution merging strategy:
#' The algorithm follows the general principles of classical multiresolution
#' segmentation:
#'
#' \enumerate{
#' \item for every region, all adjacent neighbors are evaluated;
#'
#' \item the neighboring region producing the minimum merged heterogeneity is
#' selected;
#'
#' \item a merge candidate is accepted only if both regions select each other as their best-fitting neighbor:
#'
#' \deqn{
#' A \rightarrow B
#' \quad\text{and}\quad
#' B \rightarrow A
#' }
#'
#' \item all admissible mutual pairs are sorted by increasing heterogeneity;
#'
#' \item merges are executed greedily while preventing conflicts (region can
#' participate in at most one merge per iteration);
#'
#' \item region statistics and adjacency relationships are updated
#' incrementally after each accepted merge;
#'
#' \item the process repeats until no admissible mutual merge satisfies the
#' scale criterion.
#' }
#'
#'
#' @section Merge criterion:
#' For each pair of adjacent regions, the algorithm evaluates the
#' heterogeneity of the hypothetical merged object:
#'
#' \deqn{
#' h_r =
#' SSD_1 +
#' SSD_2 +
#' \frac{W_1 W_2}{W_1 + W_2}
#' (\mu_1 - \mu_2)^2
#' }
#'
#' where:
#'
#' \itemize{
#' \item \eqn{W_1} and \eqn{W_2} are the total area weights of the two regions;
#' \item \eqn{\mu_1} and \eqn{\mu_2} are their weighted mean radiometric values;
#' \item \eqn{SSD_1} and \eqn{SSD_2} are their weighted sums of squared
#' deviations.
#' }
#'
#' This quantity corresponds to the weighted within-object heterogeneity of
#' the merged region and is evaluated incrementally without revisiting
#' individual pixels.
#'
#' Optionally, heterogeneity may be normalized by the length of the shared
#' boundary between adjacent regions:
#'
#' \deqn{
#' h =
#' \frac{h_r}
#' {B^{\alpha}}
#' }
#'
#' where \eqn{B} is the shared boundary length and \eqn{\alpha} is controlled
#' by `boundary_weight`.
#'
#' Larger values of `boundary_weight` favor merges between regions sharing
#' longer borders, whereas smaller values place greater emphasis on
#' radiometric similarity.
#'
#' Special cases are:
#'
#' \itemize{
#' \item `boundary_weight = 0`: purely radiometric merging;
#' \item `boundary_weight = 1`: heterogeneity normalized by shared boundary
#' length.
#' }
#'
#' Boundary length is estimated from boundary pixels and their associated area
#' weights.
#'
#' Boundary length is computed as the sum of shared boundary lengths between
#' neighboring pixels. The contribution of each interface is proportional to the
#' characteristic linear dimension of the corresponding spherical pixel,
#' approximated by the square root of its area.
#'
#' This quantity is a shape descriptor derived from spherical area weights. It
#' should not be interpreted as an exact geodesic perimeter on the hemisphere.
#'
#' @section Scale parameter normalization:
#'
#' The raw weighted SSD values can become numerically large because
#' heterogeneity accumulates naturally with object size.
#'
#' To improve interpretability, the user-provided `scale_parameter` is
#' normalized relative to the initial segmentation.
#'
#' First, heterogeneity is computed for all adjacent region pairs in the
#' initial segmentation. The median of these values defines a characteristic
#' heterogeneity scale:
#'
#' \deqn{
#' S =
#' median(SSD_{merge}^{initial})
#' }
#'
#' The effective internal threshold becomes:
#'
#' \deqn{
#' scale_{internal}
#' =
#' scale_{parameter} \times S
#' }
#'
#' Consequently:
#'
#' \itemize{
#' \item `scale_parameter = 1` approximately corresponds to a typical initial
#' merge;
#' \item smaller values produce finer segmentations;
#' \item larger values produce coarser segmentations.
#' }
#'
#' This normalization makes the scale parameter substantially more stable
#' across:
#'
#' \itemize{
#' \item image dynamic range;
#' \item lens projection;
#' \item spatial resolution;
#' \item spherical weighting schemes.
#' }
#'
#' @param r numeric [terra::SpatRaster-class] with one layer. Pixel values used
#'   to compute object statistics.
#' @param w [terra::SpatRaster-class] generated with [area_weight_image()].
#' @param seg numeric [terra::SpatRaster-class] with one layer and values
#'   interpreted as positive integers. Segmentation map where each pixel stores
#'   a segment identifier.
#' @param scale_parameter numeric vector of length one. Maximum admissible merge
#'   cost relative to the initial segmentation scale.
#' @param boundary_weight numeric vector of length one. Controlls the influence
#'   of shared boundary length on the merge criterion. A value of `0` produces
#'   purely radiometric merging, whereas larger values increasingly favor merges
#'   between regions sharing longer borders.
#' @param ignore_label_0 logical of length one. If `TRUE`, the segment labeled
#'   `0` is ignored.
#' @param max_iterations numeric vector of length one. Maximum number of merge
#'   iterations.
#'
#' @return Numeric [terra::SpatRaster-class] of one layer and values
#'   interpreted as positive integer and of same geometry as `seg`. Number of
#'   unique values cannot be higher than in `seg`.
#'
#' @note
#' Region statistics are updated incrementally after each merge.
#'
#' Unlike naive implementations, the algorithm does not recompute region
#' statistics from raw pixels after every iteration.
#'
#' Adjacency relationships are also updated incrementally through a dynamic
#' neighborhood graph representation.
#'
#' Multiple non-conflicting merges are executed during each iteration,
#' substantially reducing the total number of global iterations required for
#' convergence compared to purely sequential single-merge strategies.
#'
#' Region membership is also updated incrementally using cached pixel indices,
#' avoiding repeated full-image relabeling operations.
#'
#' The implementation therefore combines conceptual transparency with
#' substantially improved computational efficiency relative to naive
#' bottom-up region-merging implementations.
#'
#' @seealso [adaptive_healpix_segmentation()]
#'
#' @references \insertAllCited{}
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
#' # Initial oversegmentation
#' seg0 <- healpix_segmentation(z, a, 32)
#'
#' # Spherical pixel weights
#' w <- area_weight_image(z)
#'
#' # Multiresolution merging
#' seg <- merge_segments(
#'   r = r,
#'   w = w,
#'   seg = seg0,
#'   scale_parameter = 50
#' )
#'
#' plot(seg)
#' plot(compute_mean(r, seg))
#' display_caim(r, seg = seg)
#'
#' }
#'
#' @export
merge_segments <- function(r,
                           w,
                           seg,
                           scale_parameter,
                           boundary_weight = 0.5,
                           ignore_label_0 = TRUE,
                           max_iterations = Inf) {

  .assert_single_layer(r)
  .assert_single_layer(seg)
  .assert_area_weight_image(w)
  .assert_same_geom(r, w)
  .assert_same_geom(r, seg)
  .check_vector(scale_parameter, "numeric", 1, sign = "nonnegative")
  .check_vector(boundary_weight, "numeric", 1, sign = "nonnegative")
  .check_vector(ignore_label_0, "logical", 1)
  .check_vector(max_iterations, "numeric", 1, sign = "nonnegative")

  x_r <- terra::as.matrix(r, wide = TRUE)
  x_w <- terra::as.matrix(w, wide = TRUE)
  x_seg <- terra::as.matrix(seg, wide = TRUE)

  if (ignore_label_0) {
    x_seg[x_seg == 0] <- NA_integer_
  }

  adj <- .compute_adjacency(
    seg,
    w,
    directed = FALSE,
    ignore_label_0 = ignore_label_0,
    return_border_length = TRUE
  )

  if (!is.finite(max_iterations)) {
    max_iterations <- .Machine$integer.max
  } else {
    max_iterations <- as.integer(max_iterations)
  }

  result <- merge_segments_cpp(
    x_r = x_r,
    x_w = x_w,
    x_seg = x_seg,
    adj = adj,
    scale_parameter = scale_parameter,
    boundary_weight = boundary_weight,
    max_iterations = max_iterations
  )

  seg_out <- terra::deepcopy(seg)

  seg_out[] <- as.vector(t(result))
  seg_out[is.na(seg_out)] <- 0

  seg_out
}



