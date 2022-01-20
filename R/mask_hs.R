#' Mask hemisphere
#'
#' Given a zenith or azimuth image and angle restrictions, it produces a mask.
#'
#' @param r \linkS4class{RasterLayer}. The result of a call to
#'   \code{\link{zenith_image}} or \code{\link{azimuth_image}}.
#' @param from,to angle in degrees, inclusive limits.
#'
#' @export
#' @family Segmentation functions
#' @seealso \code{\link{masking}}
#'
#' @return An object of class \linkS4class{RasterLayer} with values \code{0} and
#'   \code{1}.
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
#' # if you want 15 degress at each side of 0
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
  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(from) == "numeric")
  stopifnot(class(to) == "numeric")
  stopifnot(length(from) == 1)
  stopifnot(length(to) == 1)

  m <- is.na(r)
  r[is.na(r)] <- 0
  r[r >= from & r <= to] <- NA
  r <- is.na(r)
  r[m] <- 0
  r
}
