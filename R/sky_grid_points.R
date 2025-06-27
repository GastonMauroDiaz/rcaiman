#' Cell central points of a sky grid
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams zenith_image
#'
#' @returns data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' }
sky_grid_points <- function(z, a, angle_width, lens_coef = lens()) {
  az <- expand.grid(seq(0 + angle_width/2, 360 - angle_width/2, angle_width),
                    seq(0 + angle_width/2, 90 - angle_width/2, angle_width))
  rr <- calc_relative_radius(az[,2], lens_coef)

  pol <- data.frame(theta = az[,1] * pi/180 + pi/2,
                    z = rr * 90 * pi/180)
  cart <- pracma::pol2cart(as.matrix(pol))
  p <- terra::vect(matrix(cart, ncol = 2))
  terra::crs(p) <- terra::crs(z)

  z <- terra::deepcopy(z)
  terra::ext(z) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
  i <- terra::rasterize(p, z)
  i <- !is.na(i)
  i <- cells(i)[i[]]
  sky_points <- terra::rowColFromCell(z, i)
  sky_points <- as.data.frame(sky_points)
  names(sky_points) <- c("row", "col")
  sky_points
}
