#' Row and col numbers from zenith and azimuth angles
#'
#' @inheritParams ootb_mblt
#' @inheritParams zenith_image
#' @param za Numeric vector of length two. Zenith and azimuth angles in degrees.
#'
#' @export
#'
#' @family HSP Functions
#' @return Numeric vector of length two.
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' row_col_from_zenith_azimuth(z, c(45, 270), lens())
row_col_from_zenith_azimuth <- function(z, za, lens_coef) {
  .is_single_layer_raster(z, "z")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(is.numeric(lens_coef))
  stopifnot(is.numeric(za))
  stopifnot(length(za) == 2)
  az <- rev(za)
  rr <- calc_relative_radius(az[2], lens_coef)
  pol <- data.frame(theta = az[1] * pi/180 + pi/2,
                    z = rr * 90 * pi/180,
                    z = 0)
  cart <- pracma::pol2cart(as.matrix(pol))
  p <- terra::vect(matrix(cart[1:2], ncol = 2))
  terra::crs(p) <- terra::crs(z)
  e <- terra::ext(z)
  terra::ext(z) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
  ir <- terra::rasterize(p, z)
  l <- list(za, terra::cells(z, terra::ext(p)) %>%
              terra::rowColFromCell(z, .) %>% as.numeric())
  names(l) <- c("zenith_azimuth", "row_col")
  l
}
