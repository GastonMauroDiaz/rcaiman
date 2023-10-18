#' Correct vignetting effect
#'
#' @inheritParams fisheye_to_equidistant
#' @param lens_coef_v Numeric vector. Coefficients of a vignetting
#'   function (\eqn{f_v}) of the type \eqn{f_v = 1 + a \cdot \theta +
#'   b \cdot \theta^2 + ... + m \cdot \theta^n}, where \eqn{\theta} is the
#'   zenith angle, \eqn{a, b, c} and \eqn{m} are the coefficients. The maximum
#'   polynomial degree supported is fifth. See [extract_radiometry()] for
#'   additional details.
#'
#' @return The argument `r` but with corrected values.
#' @export
#'
#' @family Tool Functions
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' r <- gbc(caim$Blue)
#' r
#' r <- correct_vignetting(r, z, c(0.0638, -0.101))
#' r
correct_vignetting <- function(r, z, lens_coef_v) {
  # only to avoid note from check, code is OK without this line.
  a <- b <- d <- e <- NA

  .fv <- function(theta, lens_coef_v) {
    x <- lens_coef_v[1:5]
    x[is.na(x)] <- 0
    for (i in 1:5) assign(letters[i], x[i])
    1 + a * theta + b * theta^2 + c * theta^3 + d * theta^4 + e * theta^5
  }
  r <- r / .fv(z * pi / 180, lens_coef_v)
  r[is.na(z)] <- 0
  r
}
