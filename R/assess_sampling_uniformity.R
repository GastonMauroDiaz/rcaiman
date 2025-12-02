#' Assess sampling uniformity
#'
#' @description
#' Quantify how uniformly a set of sky sampling points covers the hemispherical
#' domain with the Kullback–Leibler divergence.
#'
#' @details
#' The method projects the input points onto an equal–area grid
#' (`equalarea_seg`), computes empirical cell frequencies, and
#' evaluates the resulting discrete distribution against a uniform target.
#'
#' Let `N` be the number of valid cells and `c_i` the count in cell `i`.
#' The empirical distribution is
#' \deqn{p_i = \frac{c_i}{\sum_{j=1}^N c_j}.}
#'
#' The uniform reference is
#' \deqn{u_i = \frac{1}{N}.}
#'
#' Kullback–Leibler divergence from uniformity is
#' \deqn{\mathrm{KL}(p\|u) = \sum_i p_i \log(p_i N).}
#'
#'
#' @param equalarea_seg single-layer [terra::SpatRaster-class] with integer values. Sky
#'   segmentation map produced by [equalarea_segmentation()].
#'
#' @inheritParams extract_dn
#'
#' @return Numeric vector of length one. Kullback–Leibler divergence from uniformity.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' seg <- equalarea_segmentation(z, a, 20)
#' bin <- binarize_by_region(r, seg, method = "thr_twocorner")
#' bin <- bin & select_sky_region(z, 0, 80)
#'
#' g <- skygrid_segmentation(z, a, 5, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_black = 3)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' n_cells <- extract_dn(bin, fibonacci_points(z, a, 1), use_window = FALSE)[,3] %>% sum()
#' seg <- equalarea_segmentation(z, a, n_cells)
#' kl <- assess_sampling_uniformity(sky_points, seg)
#'
#' seg <- equalarea_segmentation(z, a, 200)
#' kl <- lapply(1:30, function(i) {
#'   sky_points <-  rem_nearby_points(sky_points, NULL, z, a, i, space = "spherical")
#'   assess_sampling_uniformity(sky_points, seg)
#' })
#' plot(unlist(kl))
#' }
assess_sampling_uniformity <- function(sky_points, equalarea_seg) {

  .assert_single_layer(equalarea_seg)
  if (names(equalarea_seg) != "Sky segments of equal area") {
    stop("`equalarea_seg` must be the output of `equalarea_segmentation()")
  }

  cell_ids <- extract_dn(equalarea_seg, sky_points, use_window = FALSE)[,3]
  if (any(is.na(cell_ids))) stop("Sky point outside the segmentation.")

  counts <- table(cell_ids)
  N <- terra::minmax(equalarea_seg)[2]

  # empirical distribution, p_i
  p <- numeric(length = N)
  names(p) <- 1:N
  p[names(counts)] <- as.numeric(counts)
  p <- p / sum(p)

  # Since the uniforme distribution is u_i = 1 / N, KL is:
  sum(p[p > 0] * log(p[p > 0] * N))

  # list(kl = kl, coverage = sum(p > 0) / N)
}
