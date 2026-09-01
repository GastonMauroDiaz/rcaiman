#' Assess sampling uniformity
#'
#' @description
#' Quantify how a set of sky sampling points covers an equal-area hemispherical
#' segmentation using (i) Kullback–Leibler divergence from a uniform
#' distribution and (ii) the fraction of cells receiving at
#' least one point, hereafter referred to as coverage.
#'
#' @details
#' Points are projected onto an equal-area segmentation (`equal_area_seg`), and
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
#' @param equal_area_seg [terra::SpatRaster-class] generated with
#'   [equal_area_segmentation()].
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
#' z <- zenith_image(1000, equidistant_coef)
#' a <- azimuth_image(z)
#'
#' row <- col <- sample(1:1000, 100)
#' pts <- data.frame(row, col)
#' i <- !is.na(extract_dn(z, pts)[,3])
#' pts <- add_directions(pts[i,], z, a)
#'
#' seg <- equal_area_segmentation(z, a, 200)
#'
#' assess_sampling_uniformity(pts, seg)
#' }
assess_sampling_uniformity <- function(sampling_points, equal_area_seg) {

  .check_columns(sampling_points, required_columns = c("row", "col"))
  .assert_equal_area_seg(equal_area_seg)

  cell_ids <- extract_dn(equal_area_seg, sampling_points, use_window = FALSE)[,3]
  if (any(is.na(cell_ids))) stop("Sky point outside the segmentation.")

  counts <- table(cell_ids)
  N <- .get_max(equal_area_seg)

  # empirical distribution, p_i
  p <- numeric(length = N)
  names(p) <- 1:N
  p[names(counts)] <- as.numeric(counts)
  p <- p / sum(p)

  # Since the uniforme distribution is u_i = 1 / N, KL is:
  kl <- sum(p[p > 0] * log(p[p > 0] * N))

  list(kl = kl, coverage = sum(p > 0) / N)
}
