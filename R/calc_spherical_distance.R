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
#' Directions are compared element by element. Within each direction set,
#' `zen` and `azi` must have the same length. The two sets must either have the
#' same length or one set must contain exactly one direction, which is then
#' compared with every direction in the other set. Other forms of vector
#' recycling are rejected. This function does not calculate all pairwise
#' combinations between two direction sets.
#'
#' Missing values are propagated element by element. Empty paired inputs return
#' `numeric(0)`.
#'
#' @param zen_1,azi_1 numeric vectors of equal length. Zenith and azimuth angles
#'   in radians defining the first direction set.
#' @param zen_2,azi_2 numeric vectors of equal length. Zenith and azimuth angles
#'   in radians defining the second direction set.
#'
#' @returns Numeric vector containing one spherical distance in radians per
#'   paired comparison.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' zen_1 <- rnorm(10, 45, 20) * pi/180
#' azi_1 <- rnorm(10, 180, 90) * pi/180
#' calc_spherical_distance(zen_1, azi_1, 0, 0)
#'
#' zen_2 <- zen_1 + 0.01
#' azi_2 <- azi_1 - 0.02
#' calc_spherical_distance(zen_1, azi_1, zen_2, azi_2)
calc_spherical_distance  <- function(zen_1, azi_1, zen_2, azi_2) {
  .check_vector(
    zen_1,
    "numeric",
    sign = "any",
    allow_na = TRUE,
    name = "zen_1"
  )
  .check_vector(
    azi_1,
    "numeric",
    sign = "any",
    allow_na = TRUE,
    name = "azi_1"
  )
  .check_vector(
    zen_2,
    "numeric",
    sign = "any",
    allow_na = TRUE,
    name = "zen_2"
  )
  .check_vector(
    azi_2,
    "numeric",
    sign = "any",
    allow_na = TRUE,
    name = "azi_2"
  )

  if (length(zen_1) != length(azi_1)) {
    stop("`zen_1` and `azi_1` must have the same length.", call. = FALSE)
  }
  if (length(zen_2) != length(azi_2)) {
    stop("`zen_2` and `azi_2` must have the same length.", call. = FALSE)
  }

  length_1 <- length(zen_1)
  length_2 <- length(zen_2)
  compatible_lengths <- length_1 == length_2 ||
    length_1 == 1L ||
    length_2 == 1L
  if (!compatible_lengths) {
    stop(
      "The two direction sets must have the same length, or one set must ",
      "contain exactly one direction.",
      call. = FALSE
    )
  }

  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  acos(pmax(pmin(cos(zen_1) * cos(zen_2) +
                   sin(zen_1) * sin(zen_2) * cos(abs(azi_2 - azi_1)), 1), -1)) %>% unname()
}
calc_spherical_distance <- compiler::cmpfun(calc_spherical_distance)
