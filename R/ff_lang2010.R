#' Generate flat-field using the Lang (2010) empirical aperture model
#'
#' @description
#' Generate a flat-field image using the aperture-dependent empirical
#' vignetting model described by \insertCite{Lang2010;textual}{rcaiman}.
#'
#' @details
#' This function implements Eq. (5) from \insertCite{Lang2010;textual}{rcaiman}.
#' The relative pixel response with respect to the image center is:
#'
#' \deqn{
#' B = 0.622 - \frac{1.504 \theta_v}{A^{4.553}} +
#' 0.0539 \ln\left((22.87 - A)(91.0 - \theta_v)\right)
#' }
#'
#' where \eqn{\theta_v} is the zenith angle in degrees (0–90 deg) and
#' \eqn{A} is the aperture value used during image acquisition.
#'
#' @param aperture numeric scalar. Aperture value used during image acquisition.
#'
#' @inheritParams ff_simple_polynomial
#'
#' @return Numeric [terra::SpatRaster-class] with the same geometry as `z`.
#'   Values represent relative directional sensitivity.
#'
#' @references \insertAllCited{}
#'
#' @export
ff_lang2010 <- function(z, aperture) {

  .assert_single_layer(z)
  .check_vector(aperture, "numeric", 1, sign = "positive")

  if (aperture > 22.87)
    stop("`aperture` value not supported.", call. = FALSE)

  ff <- 0.622 -
    (1.504 * z) / aperture^4.553 +
    0.0539 * log((22.87 - aperture) * (91.0 - z))

  vals <- terra::values(ff, mat = FALSE)

  if (any(vals[!is.na(vals)] <= 0)) {
    stop("`flat_field` must contain strictly positive values.",
         call. = FALSE)
  }
  names(ff) <- "flat_field"
  ff
}
