#' Calculate spherical distance
#'
#' @description
#' Compute the angular distance between directions on the unit sphere from
#' zenith and azimuth angles.
#'
#' @details
#' Use the spherical law of cosines to compute the angle \eqn{\gamma} between two
#' directions \eqn{(z_1, a_1)} and \eqn{(z_2, a_2)}:
#'
#' \deqn{
#' \gamma = \arccos\left(
#' \cos z_1 \cos z_2 +
#' \sin z_1 \sin z_2 \cos(a_2 - a_1)
#' \right)
#' }
#'
#' This formulation is mathematically equivalent to the dot-product-based
#' approach used internally in other functions (see [nn_spherical()]), where
#' directions are first mapped to unit Cartesian vectors and distances are
#' computed as the arc-cosine of the dot-product.
#'
#' This function is kept for clarity and transparency, as it directly expresses
#' the spherical geometry. Internally, more efficient implementations may be
#' used without changing the underlying distance definition.
#'
#' @param zen_1 numeric vector. Zenithal angle in radians.
#' @param azi_1 numeric vector. Azimuthal angle in radians.
#' @param zen_2 numeric vector of length one. Zenithal angle in radians.
#' @param azi_2 numeric vector of length one. Azimuthal angle in radians.
#'
#' @returns Numeric vector of the same length as `zen_1` and `azi_1`, containing the
#'   spherical distance (radians) to (`zen_2`, `azi_2`).
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' zen_1 <- rnorm(10, 45, 20) * pi/180
#' azi_1 <- rnorm(10, 180, 90) * pi/180
#' calc_spherical_distance(zen_1, azi_1, 0, 0)
calc_spherical_distance  <- function(zen_1, azi_1, zen_2, azi_2) {
  .check_vector(zen_1, "numeric", sign = "any", allow_na = TRUE)
  .check_vector(azi_1, "numeric", sign = "any", allow_na = TRUE)
  .check_vector(zen_2, "numeric", 1, sign = "any")
  .check_vector(azi_2, "numeric", 1, sign = "any")

  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  acos(pmax(pmin(cos(zen_1) * cos(zen_2) +
                   sin(zen_1) * sin(zen_2) * cos(abs(azi_2 - azi_1)), 1), -1)) %>% unname()
}
calc_spherical_distance <- compiler::cmpfun(calc_spherical_distance)
