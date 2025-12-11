#' Optimize the minimum distance between sky sampling points
#'
#' @description
#' Evaluate a sequence of candidate minimum distances and select the one that
#' minimizes the Kullback–Leibler (KL) divergence from a uniform distribution on
#' an equal-area hemispherical segmentation.
#'
#' @details
#' For each value in `min_dist_seq`, the input point pattern (`sky_points`) is
#' thinned using [rem_nearby_points()]. The resulting points are assigned to the
#' cells of `equalarea_seg`, and their empirical distribution is compared with
#' the uniform reference through KL divergence, as computed by
#' [assess_sampling_uniformity()]. As a results, it improves spatial regularity
#' by selecting an optimal minimum distance.
#'
#' @param min_dist_seq numeric vector. Candidate minimum distances (deg).
#'
#' @inheritParams extract_dn
#' @inheritParams rem_nearby_points
#' @inheritParams assess_sampling_uniformity
#' @inheritParams skygrid_segmentation
#'
#' @return List with:
#' \describe{
#'   \item{`sky_points`}{Thinned point set using the optimal distance.}
#'   \item{`min_dist`}{Value in `min_dist_seq` that minimized KL divergence.}
#'   \item{`kl_values`}{KL divergence for each candidate distance.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' bin <- binarize_by_region(r, ring_segmentation(z, 15), "thr_isodata") &
#'   select_sky_region(z, 0, 88)
#'
#' seg <- skygrid_segmentation(z, a, 10)
#' sky_points <- sample_sky_points(r, bin, seg,
#'                                 dist_to_black = 3)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' l <- optimize_sampling_uniformity(sky_points, NULL, z, a, equalarea_segmentation(z, a, 100), 0:15)
#' plot(l$kl_values)
#' l$min_dist
#'
#' plot(bin)
#' points(l$sky_points$col, nrow(caim) - l$sky_points$row, col = 2, pch = 10)
#' }
optimize_sampling_uniformity <- function(sky_points,
                                         r,
                                         z,
                                         a,
                                         equalarea_seg,
                                         min_dist_seq) {

  .check_vector(min_dist_seq, "numeric", sign = "nonnegative")
  #more basic checks are handled by the functions called below

  kl_values <- lapply(min_dist_seq, function(i) {
    sky_points <-  rem_nearby_points(sky_points, r, z, a, i, space = "spherical")
    assess_sampling_uniformity(sky_points, equalarea_seg)$kl
  })

  best <- min_dist_seq[which.min(kl_values)]
  sky_points <- rem_nearby_points(sky_points, r, z, a, min_dist = best, space =  "spherical")

  list(sky_points = sky_points,
       min_dist   = best,
       kl_values  = unlist(kl_values))
}
