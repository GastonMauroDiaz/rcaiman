#' Correct vignetting using legacy calibration model
#'
#' @description
#' Apply vignetting correction using the aperture-dependent empirical model
#' described by \insertCite{Lang2010;textual}{rcaiman}.
#'
#' @details
#' This function implements the aperture-dependent vignetting model reported in
#' Eq. (5) of \insertCite{Lang2010;textual}{rcaiman}. The relative pixel value
#' with respect to the image center is given by:
#'
#' \deqn{
#' B = 0.622 - \frac{1.504 \theta_v}{A^{4.553}} +
#' 0.0539 \ln\left((22.87 - A)(91.0 - \theta_v)\right)
#' }
#'
#' where \eqn{\theta_v} is the zenith angle in degrees (0–90 deg) and
#' \eqn{A} is the aperture value used during image acquisition.
#'
#' This function is provided to support workflows based on that historical
#' calibration. It does not use the polynomial schemes implemented in
#' [correct_vignetting()].
#'
#' @inheritParams fisheye_to_equidistant
#' @param aperture numeric vector of length one. Aperture value used during
#'   image acquisition.
#'
#' @return Numeric [terra::SpatRaster-class] with the same geometry as `r`,
#'   with pixel values divided by \eqn{B} to correct for vignetting.
#'
#' @references \insertAllCited{}
#'
#' @export
correct_vignetting_legacy <- function(r, z, aperture) {

  .check_r_z_a_m(r, z, r_type = "any")
  .check_vector(aperture, "numeric", 1, sign = "positive")

  if (aperture > 22.87) {
    stop("`aperture` value not supported.")
  }

  B <- 0.622 -
    (1.504 * z) / aperture^4.553 +
    0.0539 * log((22.87 - aperture) * (91.0 - z))

  r / B
}

