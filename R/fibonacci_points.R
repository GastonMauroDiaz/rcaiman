#' Map Fibonacci points to raster coordinates
#'
#' Return image row and column indices for points distributed on the upper
#' hemisphere using a spherical Fibonacci lattice with approximately constant
#' angular spacing.
#'
#' @param spacing numeric vector of length one. Angular separation (deg) between
#'
#' @inheritParams skygrid_segmentation
#'
#' @return `data.frame` with integer columns `row` and `col`.
#'
#' @seealso [skygrid_centers()]
#'
#' @export
#'
#' @examples
#' z <- zenith_image(100, lens())
#' a <- azimuth_image(z)
#' sampling_points <- fibonacci_points(z, a, 30)
#' \dontrun{
#' display_caim(is.na(z), sky_points = sky_points)
#' }
fibonacci_points <- function(z, a, spacing) {
  .check_r_z_a_m(NULL, z, a)
  .check_vector(spacing, "numeric", 1, sign = "positive")

  # --- internal helper: fibonacci hemisphere generation ---
  fibonacci_hemisphere <- function(spacing) {
    theta <- spacing * pi / 180
    d <- 2 * sin(theta / 2)
    N_full <- round(4 * pi / (d^2))
    N <- round(N_full / 2)
    phi <- (1 + sqrt(5)) / 2

    pts <- matrix(NA_real_, nrow = N, ncol = 3)

    for (k in seq_len(N)) {
      zc <- 1 - (k - 0.5) / N
      r <- sqrt(1 - zc * zc)
      az <- 2 * pi * ((k - 1) / phi %% 1)
      x <- r * cos(az)
      y <- r * sin(az)
      pts[k, ] <- c(x, y, zc)
    }
    pts
  }

  # --- convert Cartesian to spherical angles ---
  cartesian_to_angles <- function(xyz) {
    x <- xyz[,1]; y <- xyz[,2]; zc <- xyz[,3]
    zenith <- acos(zc)
    azimuth <- atan2(y, x)
    azimuth[azimuth < 0] <- azimuth[azimuth < 0] + 2*pi
    cbind(zenith = zenith, azimuth = azimuth)
  }

  # generate Fibonacci points on hemisphere
  pts_xyz <- fibonacci_hemisphere(spacing)
  pts_ang <- cartesian_to_angles(pts_xyz)

  # convert to degrees
  zen <- pts_ang[, "zenith"] * 180/pi
  i <- zen < 89
  zen <- zen[i]
  azi <- pts_ang[i, "azimuth"] * 180/pi


  # project to relative radius using lens model
  rr <- calc_relative_radius(zen, attr(z, "lens_coef"))

  # convert to projected Cartesian coordinates used by the image
  pol <- data.frame(
    theta = azi * pi/180 + pi/2,
    z = rr * 90 * pi/180
  )
  cart <- pracma::pol2cart(as.matrix(pol))
  p <- terra::vect(matrix(cart, ncol = 2))
  terra::crs(p) <- terra::crs(z)

  # rasterize to locate each point in pixel space
  ztmp <- terra::deepcopy(z)
  terra::ext(ztmp) <- terra::ext(-pi/2, pi/2, -pi/2, pi/2)

  r <- terra::rasterize(p, ztmp)
  r <- !is.na(r)
  idx <- cells(r)[r[]]
  ij <- terra::rowColFromCell(ztmp, idx)
  ij <- as.data.frame(ij)
  names(ij) <- c("row", "col")

  ij
}
