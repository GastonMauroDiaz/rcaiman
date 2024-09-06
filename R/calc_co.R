#' Calculate canopy openness
#'
#' Calculate canopy openness
#'
#' Canopy openness calculated as in the equation from
#' \insertCite{Gonsamo2011;textual}{rcaiman}:
#'
#' \eqn{CO = \sum_{i = 1}^{N} GF(\phi_i, \theta_i) \cdot [(cos(\theta_1) -
#' cos(\theta_2))/n]},
#'
#' where \eqn{GF(\phi_i, \theta_i)} is the gap fraction of the cell \eqn{i},
#' \eqn{\theta_1} and \eqn{\theta_2} are the minimum and maximum zenith angle of
#' the cell \eqn{i}, \eqn{n} is the number of cells on the ring delimited by
#' \eqn{\theta_1} and \eqn{\theta_2}, and \eqn{N} is the total number of cells.
#'
#' When a mask is provided through the `m` argument, the equation is modified as
#' follow:
#'
#'
#' \eqn{
#' \frac{
#' CO = \sum_{i = 1}^{N} GF(\phi_i, \theta_i) \cdot [(cos(\theta_1) -
#' cos(\theta_2))/n]
#' }{ \sum_{i = 1}^{N} (cos(\theta_1) - cos(\theta_2))/n}
#' }.
#'
#' This allows the masking of any individual cell.
#'
#' @inheritParams sky_grid_segmentation
#' @param bin [SpatRaster-class]. Binarized hemispherical canopy image.
#' @inheritParams masking
#'
#' @return Numeric vector of length one.
#' @export
#'
#' @family Metrics Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- mask_hs(z, 0, 70)
#' bin <- apply_thr(caim$Blue, thr_isodata(caim$Blue[m]))
#' plot(bin)
#' calc_co(bin, z, a, m, 10)
#'
calc_co <- function(bin, z, a, m = NULL, angle_width = 10) {
  g <- sky_grid_segmentation(z, a, angle_width)
  g[!m] <- 0
  ds <- extract_feature(bin, g, return_raster = FALSE)
  ids <- .decode_label(as.numeric(names(ds)))
  mx_angles <- ids$ring_ID * angle_width * pi/180
  mn_angles <- mx_angles - angle_width * pi/180
  .n <- 360/angle_width
  if (is.null(m)) {
    sum(ds * ((cos(mn_angles) - cos(mx_angles))/.n))
  } else {
    sum(ds * ((cos(mn_angles) - cos(mx_angles))/.n)) /
      sum((cos(mn_angles) - cos(mx_angles))/.n)
      # (cos(mn_angles %>% min()) - cos(mx_angles %>% max()))
  }
}
