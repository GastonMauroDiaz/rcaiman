#' Calculate spherical distance
#'
#' Calculate angles in radians between objects in the sky, which are
#' proportional to the shortest distances between then over a conceptual sky
#' vault.
#'
#' @param z1 Numeric vector. Zenithal angle in radians.
#' @param a1 Numeric vector. Azimuthal angle in radians.
#' @param z2 Numeric vector of length one. Zenithal angle in radians.
#' @param a2 Numeric vector of length one. Azimuthal angle in radians.
#'
#' @returns Numeric vector.
#' @export
#'
#' @examples
#' set.seed(1)
#' z1 <- rnorm(10, 45, 20) * pi/180
#' a1 <- rnorm(10, 180, 90) * pi/180
#' calc_spherical_distance(z1, a1, 0, 0)
calc_spherical_distance  <- function(z1, a1, z2, a2) {
  stopifnot(all(is.numeric(z1),is.numeric(a1), is.numeric(z2), is.numeric(a2)))
  stopifnot(length(z1) == length(a1))
  stopifnot(all(length(z2) == 1, length(a2) == 1))

  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  acos(pmax(pmin(cos(z1) * cos(z2) +
                   sin(z1) * sin(z2) * cos(abs(a2 - a1)), 1), -1))
}
