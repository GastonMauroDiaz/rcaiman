#' Mask image
#'
#' Given angle restrictions, produce an image for masking fisheye photos.
#'
#'
#' @inheritParams sky_grid_segmentation
#' @param zlim Numeric vector of length two. Angles in degrees. Set the zenith
#'   angle range with inclusive limits.
#' @param alim Numeric vector of length two. Angles in degrees. Set the azimuth
#'   angle range with inclusive limits.
#'
#' @return \linkS4class{RasterLayer}
#' @export
#'
#' @family masking functions
#' @seealso \code{\link{write_bin}}
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#'
#' m <- mask_image(z, a, c(20, 70), c(90, 180))
#' plot(m)
#'
#' m1 <- mask_image(z, a, alim = c(90, 180))
#' plot(m1)
#'
#' m2 <- mask_image(z, zlim = c(20, 70))
#' plot(m2)
#'
#' plot(m1 & m2)
#'
#' m <- mask_image(z)
#' plot(m)
mask_image <- function(z,
                       a = NULL,
                       zlim = NULL,
                       alim = NULL) {
  no_data_area <- is.na(z)

  if (all(is.null(zlim), is.null(alim))) {
    m <- !is.na(z)
  } else {
    if (all(!is.null(zlim), !is.null(alim))) {
      stopifnot(length(zlim) == 2)
      stopifnot(length(alim) == 2)

      stopifnot(all(zlim[1] >= 0, zlim[2] <= 90))
      stopifnot(all(alim[1] >= 0, alim[2] <= 360))

      z[is.na(z)] <- 0
      a[is.na(a)] <- 0
      z[z >= zlim[1] & z <= zlim[2]] <- NA
      a[a >= alim[1] & a <= alim[2]] <- NA
      m <- is.na(z) + is.na(a)
      m[m == 2] <- NA
      m <- is.na(m)
    } else {
      if (!is.null(zlim)) {
        stopifnot(length(zlim) == 2)
        stopifnot(all(zlim[1] >= 0, zlim[2] <= 90))

        z[is.na(z)] <- 0

        z[z >= zlim[1] & z <= zlim[2]] <- NA
        m <- is.na(z)
      } else {
        stopifnot(length(alim) == 2)
        stopifnot(all(alim[1] >= 0, alim[2] <= 360))

        a[is.na(a)] <- 0
        a[a >= alim[1] & a <= alim[2]] <- NA
        m <- is.na(a)
      }
    }
  }
  # fix inclusion of the area outside the circle if zmin is 0
  m[no_data_area] <- 0
  m
}
