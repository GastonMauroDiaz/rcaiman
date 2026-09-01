#' Generate flat-field with a free-form polynomial model
#'
#' @description
#' Generate a flat-field image from zenith angle data using an unconstrained
#' polynomial vignetting model.
#'
#' @details
#' The model is defined as
#'
#' \deqn{
#' f_v(\theta) =
#' a + b\theta + c\theta^2 +
#' d\theta^3 + e\theta^4 +
#' f\theta^5 + g\theta^6
#' }
#'
#' where \eqn{\theta} is the zenith angle in radians and
#' \eqn{a,b,\dots,g} are polynomial coefficients. Polynomial degrees up to 6
#' (7 parameters including the constant term) are supported. Missing
#' higher-order terms are assumed to be zero.
#'
#' This formulation does not impose structural normalization constraints. It
#' assumes that radiometry data used to obtain the coefficient is off high
#' accuracy and completeness.
#'
#' @param parameters numeric vector of length between 1 and 7.
#'   Polynomial coefficients ordered as \eqn{(a, b, c, d, e, f, g)}.
#'
#' @inheritParams ff_simple_polynomial
#'
#' @return Numeric [terra::SpatRaster-class] with one layer and the same
#'   geometry as `z`. Values represent relative directional sensitivity.
#'
#' @export
ff_free_form_polynomial <- function(z, parameters) {

  .assert_single_layer(z)
  .check_vector(parameters, "numeric", sign = "any")

  if (length(parameters) == 0 || length(parameters) > 7)
    stop("`parameters` must have length between 1 and 7.", call. = FALSE)

  # Pad parameters up to degree 6
  x <- rep(0, 7)
  x[seq_along(parameters)] <- parameters
  x[is.na(x)] <- 0

  theta <- z * pi / 180

  ff <- x[1] +
    x[2] * theta +
    x[3] * theta^2 +
    x[4] * theta^3 +
    x[5] * theta^4 +
    x[6] * theta^5 +
    x[7] * theta^6

  # Physical validation
  vals <- terra::values(ff, mat = FALSE)

  if (any(vals[!is.na(vals)] <= 0)) {
    stop("`flat_field` must contain strictly positive values.",
         call. = FALSE)
  }

  names(ff) <- "flat_field"
  ff
}
