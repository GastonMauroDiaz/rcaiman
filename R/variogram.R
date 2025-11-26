spherical_variogram <- function(z, a, values, n_lags = 20) {
  # z, a: zenith and azimuth in radians (numeric vectors, same length)
  # values: numeric vector with one value per point
  # n_lags: number of lag bins

  n <- length(z)
  idx <- utils::combn(n, 2)

  d <- calc_spherical_distance(
    z1 = z[idx[1, ]],
    a1 = a[idx[1, ]],
    z2 = z[idx[2, ]],
    a2 = a[idx[2, ]]
  )

  # semivariance
  v  <- values[idx[1, ]] - values[idx[2, ]]
  v2 <- 0.5 * v^2

  # lag bins
  breaks <- seq(0, max(d), length.out = n_lags + 1)
  bins <- cut(d, breaks = breaks, include.lowest = TRUE)

  data.frame(
    lag   = tapply(d, bins, mean),
    gamma = tapply(v2, bins, mean),
    n     = tapply(v2, bins, length)
  )
}
