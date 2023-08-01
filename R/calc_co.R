#' Calculate canopy openness
#'
#' Calculate canopy openness
#'
#' Canopy openness calculated as in the equation from
#' \insertCite{Gonsamo2011;textual}{rcaiman}:
#'
#' \eqn{CO = \sum_{i = 1}^{N} GF(\phi_i, \theta_i) \cdot ((cos(\theta_1) -
#' cos(\theta_2))/n)}, where \eqn{GF(\phi_i, \theta_i)} is the gap fraction of
#' the cell \eqn{i}, \eqn{\theta_1} and \eqn{\theta_2} are the smallest and
#' largest zenith angle of the cell \eqn{i}, and \eqn{n} is the number of cells
#' on the ring delimited by \eqn{\theta_1} and \eqn{\theta_2}.
#'
#'
#' @inheritParams sky_grid_segmentation
#' @param bin \linkS4class{SpatRaster}. Binarized hemispherical canopy image.
#' @param m \linkS4class{SpatRaster}. Check \code{\link{mask_hs}}.
#'
#' @return Numeric vector of length one.
#' @export
#'
#' @family Tool Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' bin <- apply_thr(caim$Blue, thr_isodata(caim$Blue[m]))
#' plot(bin)
#' calc_co(bin, m, z, a, 10)
#'
calc_co <- function(bin, m, z, a, angle_width) {
  g <- sky_grid_segmentation(z, a, angle_width)
  g[!m] <- NA
  ds <- extract_feature(bin, g, return_raster = FALSE)

  ids <- .decode_label(as.numeric(names(ds)))
  rcl <- tapply(ids$sector_ID, ids$ring_ID, length)
  rcl <- data.frame(ring_id = names(rcl), no = rcl)
  .n <- ids$ring_ID
  for (i in seq_along(rcl$ring_id)) {
    .n[.n == rcl$ring_id[i]] <- rcl$no[i]
  }
  mx_angles <- ids$ring_ID * angle_width * pi/180
  mn_angles <- mx_angles - angle_width * pi/180
  sum(ds * ((cos(mn_angles) - cos(mx_angles))/.n))
}
