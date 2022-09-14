#' Zenith and azimuth angles from row and col numbers
#'
#' @inheritParams ootb_mblt
#' @param row_col Numeric vector of length two. Row and col numbers.
#' @inheritParams zenith_image
#'
#' @export
#'
#' @family HSP Functions
#' @return Numeric vector of length two.
#' @examples
#' z <- zenith_image(1000, lens_coef = lens())
#' zenith_azimuth_from_row_col(z, c(501, 750), lens())
zenith_azimuth_from_row_col <- function(z, row_col, lens_coef) {
  .is_single_layer_raster(z, "z")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(is.numeric(lens_coef))
  stopifnot(is.numeric(row_col))
  stopifnot(length(row_col) == 2)

  zz <- terra::rast(z)

  #get azimuth
  e <- terra::ext(zz)
  terra::ext(zz) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
  xy <- terra::cellFromRowCol(zz, row_col[1], row_col[2]) %>%
    terra::xyFromCell(zz, .)
  tr <- pracma::cart2pol(as.numeric(xy))
  azimuth <- tr[1] - pi/2 * 180/pi
  if (azimuth < 0) azimuth <- 360 + azimuth
  #get relative radius
  rr <- tr[2] * 180/pi / 90
  #invert
  zs <- seq(0,150, 0.1)
  rrs <- calc_relative_radius(zs, lens_coef)
  z_from_rr <- suppressWarnings(splinefun(rrs, zs))
  zenith <- z_from_rr(rr)
  l <- list(c(zenith, azimuth), row_col)
  names(l) <- c("zenith_azimuth", "row_col")
  l
}
