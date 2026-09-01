#' Build area weight image
#'
#' @description
#' Build a single-layer image with relative solid angle weights.
#'
#' @details
#' In an ideal equal-area projection, every pixel represents the same solid
#' angle of the hemisphere. Real lenses deviate from this behavior due to
#' geometric distortion. This function computes relative solid angle weights
#' analytically from the lens projection model.
#'
#' The relative solid angle weight \eqn{W} is calculated analytically as:
#' \deqn{
#' W(\theta) = \frac{\sin(\theta)}{f(\theta) \cdot f'(\theta)}
#' }
#' where \eqn{f(\theta)} is the lens projection polynomial and
#' \eqn{f'(\theta)} is its first derivative.
#'
#' At the zenith (\eqn{\theta = 0}), the analytical limit converges to:
#' \deqn{W(0) = \frac{1}{a^2}}
#' Where \eqn{a} is the first coefficient of the lens polynomial.
#'
#' @inheritParams add_directions
#'
#' @note
#' To obtain pixel solid angles in steradians, divide the resulting raster by
#' the squared radius in pixels at 90 deg:
#' `(ncol(z) / 2)^2`.
#'
#' @return A `SpatRaster` of relative solid angle weights of same dimensions as `z`.
#'
#' @export
#'
#' @examples
#' z <- zenith_image(100, equidistant_coef)
#' summary(pixel_area_weight(z))
#' z <- zenith_image(100, equal_area_coef)
#' summary(pixel_area_weight(z))
area_weight_image <- function(z) {
  .assert_zenith_image(z)
  projection_coef <- attr(z, "projection_coef")
  if (length(projection_coef) == 1) projection_coef <- c(projection_coef, 0)
  theta <- z * pi/180

  # Horner method
  n <- length(projection_coef)

  f_t  <- z
  f_t[] <- 0
  f_prime_t <- f_t
  for (i in n:1) {
    f_t <- (f_t + projection_coef[i]) * theta
  }

  projection_coef_prime <- projection_coef * seq_along(projection_coef)
  projection_coef_prime <- projection_coef_prime[-1]
  n <- length(projection_coef_prime)
  for (i in n:1) {
    f_prime_t <- (f_prime_t + projection_coef_prime[i]) * theta
  }
  f_prime_t <- f_prime_t + projection_coef[1]


  # Apply analytical calculation handling the zenith singularity robustly
  weight <- terra::ifel(theta < 1e-7,
                        1 / (projection_coef[1]^2),
                        sin(theta) / (f_t * f_prime_t))
  names(weight) <- "Area weight image"
  weight
}


