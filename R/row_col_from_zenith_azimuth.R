#' Obtain row and col numbers from zenith and azimuth angles
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams zenith_image
#' @param zenith_azimuth Numeric vector of length two. Zenith and azimuth angles in degrees.
#'
#' @note
#' Use the `lens_coef` argument to calculate coordinates below the horizon.
#'
#' @export
#'
#' @family HSP Functions
#' @return Numeric vector of length two.
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' row_col_from_zenith_azimuth(z, a, c(45, 270))
row_col_from_zenith_azimuth <- function(z, a, zenith_azimuth, lens_coef = NULL) {
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(z, a)
  stopifnot(is.numeric(zenith_azimuth))
  stopifnot(length(zenith_azimuth) == 2)

  if (is.null(lens_coef)) {
    i <- (z - zenith_azimuth[1])^2 + (a - zenith_azimuth[2])^2
    row_col <- terra::rowColFromCell(z, which.min(i[]))
    l <- list(zenith_azimuth, row_col %>% as.numeric())
  } else {
    az <- rev(zenith_azimuth)
    rr <- calc_relative_radius(az[2], lens_coef)
    pol <- data.frame(theta = az[1] * pi/180 + pi/2,
                      z = rr * 90 * pi/180,
                      z = 0)
    cart <- pracma::pol2cart(as.matrix(pol))
    p <- terra::vect(matrix(cart[1:2], ncol = 2))
    terra::crs(p) <- terra::crs(z)

    z <- terra::deepcopy(z)
    terra::ext(z) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
    ir <- terra::rasterize(p, z)
    l <- list(zenith_azimuth, terra::cells(z, terra::ext(p)) %>%
                terra::rowColFromCell(z, .) %>% as.numeric())
  }
  names(l) <- c("zenith_azimuth", "row_col")
  l
}
