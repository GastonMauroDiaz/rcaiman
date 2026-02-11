#' Correct vignetting effect
#'
#' Apply a vignetting correction to an image using a polynomial model.
#'
#' Vignetting is the gradual reduction of image brightness toward the periphery.
#' This function applies a parametric correction tailored to a specific
#' hemispherical camera system that varies with zenith angle at the pixel level.
#'
#' The selected calibration scheme defines the modeling conventions used by
#' `rcaiman` to represent and correct vignetting effects. Different schemes
#' correspond to different assumptions and degrees of freedom in the underlying
#' polynomial model.
#'
#' \describe{
#'   \item{simple}{
#'   A constrained polynomial scheme in which the vignetting function is
#'   normalized such that \eqn{f_v(0) = 1}. The effective model is
#'   \eqn{f_v(\theta) = 1 + a\theta + b\theta^2 + \dots + m\theta^n}, where
#'   \eqn{\theta} is the zenith angle (in radians) and \eqn{a,b,\dots,m} are the
#'   polynomial coefficients. Polynomial degrees up to 6 are supported.
#'   This scheme follows the approach described in
#'   \insertCite{Diaz2024;textual}{rcaiman}. See [extract_radiometry()] for guidance
#'   on estimating these coefficients.}
#'
#'   \item{free_form}{
#'   A flexible polynomial scheme in which all model parameters, including the
#'   constant term, are estimated from data. The effective model is
#'   \eqn{f_v(\theta) = a + b\theta + c\theta^2 + \dots + m\theta^n}, where
#'   \eqn{\theta} is the zenith angle (in radians) and \eqn{a,b,\dots,m} are free
#'   parameters. Up to 7 parameters are supported.
#'   This scheme assumes that the input data reliably capture the true shape of
#'   the vignetting function, and therefore applies no structural constraints on
#'   the model parameters. It is typically used with data acquired using a
#'   photometric sphere, following the methodology described in
#'   \insertCite{Lang2010;textual}{rcaiman}.}
#' }
#'
#' @inheritParams fisheye_to_equidistant
#'
#' @param lens_coef_v numeric vector. Coefficients of the vignetting function.
#'   See *Details*.
#' @param scheme character vector of length one. Calibration scheme defining
#'   the modeling assumptions used by `rcaiman` for radiometric correction.
#'   The selected scheme determines the effective form of the model and the
#'   constraints applied to its parameters. Supported values are
#'   `"simple"` and `"free_form"`.
#'
#' @param model_type character vector of lenght one. Only `"polynomial"` is currently
#'   supported.
#'
#' @return [terra::SpatRaster-class] with the same content as `r` but with
#'   pixel values adjusted to correct for vignetting, preserving all other
#'   properties (layers, names, extent, and CRS).
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/APC_0836.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(2132, lens("Olloclip"))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)
#'
#' caim <- expand_noncircular(caim, z, zenith_colrow)
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#'
#' bin <- binarize_with_thr(caim$Blue, thr_isodata(caim$Blue[m]))
#' display_caim(caim$Blue, bin)
#'
#' caim <- invert_gamma_correction(caim, 2.2)
#' caim <- correct_vignetting(caim, z, c(-0.0546, -0.561, 0.22)) %>%
#'         normalize_minmax()

# The lens_coef_v values are from doi:10.1016/j.agrformet.2024.110020
#' }
correct_vignetting <- function(r, z, lens_coef_v, method = "simple", model_type = "polynomial") {

  .check_r_z_a_m(r, z, r_type = "any")
  .check_vector(lens_coef_v, "numeric", sign = "any")
  .assert_choice(method, c("simple", "photometric_sphere"))
  .assert_choice(model_type, "polynomial")

  if (method == "simple") {
    # only to avoid note from check, code is OK without this line.
    a <- b <- d <- e <- f <- NA

    .fv <- function(theta, lens_coef_v) {
      x <- lens_coef_v[1:6]
      x[is.na(x)] <- 0
      for (i in 1:6) assign(letters[i], x[i])
      1 + a * theta + b * theta^2 + c * theta^3 +
        d * theta^4 + e * theta^5 + f * theta^6
    }
  } else {
    a <- b <- d <- e <- f <- g <- NA

    .fv <- function(theta, lens_coef_v) {
      x <- lens_coef_v[1:7]
      x[is.na(x)] <- 0
      for (i in 1:7) assign(letters[i], x[i])
      a + b * theta + c * theta^2 + d * theta^3 +
        e * theta^4 + f * theta^5 + g * theta^6
    }

  }

  r <- r / .fv(z * pi / 180, lens_coef_v)
  r[is.na(z)] <- 0
  r
}
