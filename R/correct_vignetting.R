#' Correct vignetting effect
#'
#' @inheritParams fisheye_to_equidistant
#' @param lens_coef_v Numeric vector. Coefficients of a vignetting
#'   function (\eqn{f_v}) of the type \eqn{f_v = 1 + a \cdot \theta +
#'   b \cdot \theta^2 + ... + m \cdot \theta^n}, where \eqn{\theta} is the
#'   zenith angle, \eqn{a, b, c} and \eqn{m} are the coefficients. The maximum
#'   polynomial degree supported is sixth. See [extract_radiometry()] for
#'   additional details.
#'
#' @return The argument `r` but with corrected values.
#' @export
#'
#' @family Tool Functions
#' @examples
#' \dontrun{
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(2132/2,  c(0.7836, 0.1512, -0.1558))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)/2
#'
#' caim <- expand_noncircular(caim, z, zenith_colrow)
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#'
#' bin <- apply_thr(caim$Blue, thr_isodata(caim$Blue[m]))
#'
#' display_caim(caim$Blue, bin)
#'
#' caim <- gbc(caim, 2.2)
#' caim <- correct_vignetting(caim, z, c(-0.0546, -0.561, 0.22)) %>%
#'                                                      normalize_minmax()
#' # The lens_coef_v are here coef 10.1016/j.agrformet.2024.110020
#' }
correct_vignetting <- function(r, z, lens_coef_v) {
  # only to avoid note from check, code is OK without this line.
  a <- b <- d <- e <- f <- NA

  .fv <- function(theta, lens_coef_v) {
    x <- lens_coef_v[1:6]
    x[is.na(x)] <- 0
    for (i in 1:6) assign(letters[i], x[i])
    1 + a * theta + b * theta^2 + c * theta^3 +
      d * theta^4 + e * theta^5 + f * theta^6
  }
  r <- r / .fv(z * pi / 180, lens_coef_v)
  r[is.na(z)] <- 0
  r
}
