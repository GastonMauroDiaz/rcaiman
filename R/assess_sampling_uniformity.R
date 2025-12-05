#' Assess sampling uniformity
#'
#' @description
#' Quantify how a set of sky sampling points covers an equal-area hemispherical
#' segmentation using (i) Kullback–Leibler divergence from a uniform
#' distribution and (ii) coverage, defined as the fraction of cells receiving at
#' least one point.
#'
#' @details
#' Points are projected onto an equal-area segmentation (`equalarea_seg`), and
#' cell-wise counts are computed.
#'
#' Let be \eqn{N_{\text{points}}} the total number of points, \eqn{N_{\text{cells}}}
#' the total number of cells, and \eqn{c_i} the point count in cell \eqn{i}
#' (\eqn{i = 1,\dots,N_{\text{cells}}}). The empirical distribution for cell
#' \eqn{i} is:
#'
#'   \deqn{p_i = \frac{c_i}{N_{\text{points}}}, \quad \text{with} \quad
#'   \sum_{i=1}^{N_{\text{cells}}} p_i = 1.}
#'
#' The uniform reference is \deqn{u_i = 1 / N_{\text{cells}}.}
#'
#' Kullback–Leibler divergence is computed excluding zero-probability terms:
#'   \deqn{\mathrm{KL}(p\|u) = \sum_{i: p_i > 0}^{N_{\text{cells}}} p_i \cdot \log(p_i \cdot N_{\text{cells}}).}
#'
#' Coverage is defined as:
#'   \deqn{\text{coverage} = \frac{\text{Count}(\{i : c_i > 0\})}{N_{\text{cells}}}.}
#'
#' @param equalarea_seg [terra::SpatRaster-class] generated with [equalarea_segmentation()].
#'
#' @inheritParams extract_dn
#'
#' @return List with two components:
#' \describe{
#'   \item{`kl`}{numeric value. Kullback–Leibler divergence from uniformity.}
#'   \item{`coverage`}{numeric value in \eqn{[0,1]}. Fraction of cells with at
#'     least one point.}
#' }
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
#' g <- skygrid_segmentation(z, a, 5)
#' sky_points <- sample_sky_points(r, bin, g, dist_to_black = 3)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' seg <- equalarea_segmentation(z, a, 200)
#' kl <- lapply(1:10, function(i) {
#'   sky_points <-  rem_nearby_points(sky_points, NULL, z, a, i, space = "spherical")
#'   assess_sampling_uniformity(sky_points, seg)$kl
#' })
#' plot(unlist(kl), type = "l")
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
  kl <- sum(p[p > 0] * log(p[p > 0] * N))

  list(kl = kl, coverage = sum(p > 0) / N)
}
