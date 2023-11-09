#' Mask hemisphere
#'
#' Given a zenith or azimuth image and angle restrictions, this function
#' produces a mask.
#'
#' @param r [SpatRaster-class] built with [zenith_image()] or
#'   [azimuth_image()].
#' @param from,to angle in degrees, inclusive limits.
#'
#' @export
#' @family Segmentation Functions
#' @seealso [masking()]
#'
#' @return An object of class [SpatRaster-class] with values `0` and
#'   `1`.
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' m1 <- mask_hs(z, 20, 70)
#' plot(m1)
#' m2 <- mask_hs(a, 330,360)
#' plot(m2)
#' plot(m1 & m2)
#' plot(m1 | m2)
#'
#' # 15 degrees at each side of 0
#' m1 <- mask_hs(a, 0, 15)
#' m2 <- mask_hs(a, 345, 360)
#' plot(m1 | m2)
#'
#' # better use this
#' plot(!is.na(z))
#' # instead of this
#' plot(mask_hs(z, 0, 90))
#' }
mask_hs <- function(r, from, to) {
  .is_single_layer_raster(r, "r")
  stopifnot(class(from) == "numeric")
  stopifnot(class(to) == "numeric")
  stopifnot(length(from) == 1)
  stopifnot(length(to) == 1)

  m <- is.na(r)
  r[is.na(r)] <- 0
  r[r >= from & r <= to] <- NA
  r <- is.na(r)
  r[m] <- 0
  as.logical(r)
}
