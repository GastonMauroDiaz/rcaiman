#' Select a portion of the sky vault
#'
#' Given a zenith or azimuth image and angle restrictions, this function
#' produces select a region of the sky vault and produce a binary raster mask.
#'
#' @param r [SpatRaster-class] built with [zenith_image()] or
#'   [azimuth_image()].
#' @param from,to angle in degrees, inclusive limits.
#'
#' @export
#'
#' @seealso [masking()]
#'
#' @return An object of class [SpatRaster-class] with values `0` and
#'   `1`.
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' m1 <- select_sky_vault_region(z, 20, 70)
#' plot(m1)
#' m2 <- select_sky_vault_region(a, 330,360)
#' plot(m2)
#' plot(m1 & m2)
#' plot(m1 | m2)
#'
#' # 15 degrees at each side of 0
#' m1 <- select_sky_vault_region(a, 0, 15)
#' m2 <- select_sky_vault_region(a, 345, 360)
#' plot(m1 | m2)
#'
#' # You can use this
#' plot(!is.na(z))
#' # instead of this
#' plot(select_sky_vault_region(z, 0, 90))
#' }
select_sky_vault_region <- function(r, from, to) {
  .is_single_layer_raster(r, "r")
  stopifnot(class(from) == "numeric")
  stopifnot(class(to) == "numeric")
  stopifnot(length(from) == 1)
  stopifnot(length(to) == 1)

  if (to >= max(r[], na.rm = TRUE)) {
    c1 <- !is.na(r)
  } else {
    c1 <- !apply_thr(r, to)
  }

  c1 & apply_thr(r, from)

  # m <- is.na(r)
  # r[is.na(r)] <- 0
  # r[r >= from & r <= to] <- NA
  # r <- is.na(r)
  # r[m] <- 0
  # as.logical(r)
}
