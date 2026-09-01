#' Generate flat-field with a simple polynomial model
#'
#' @description
#' Generate a flat-field image from zenith angle data using a constrained
#' polynomial vignetting model.
#'
#' @details
#' The model is normalized such that \eqn{f_v(0) = 1}. It is defined as
#'
#' \deqn{
#' f_v(\theta) =
#' 1 + a\theta + b\theta^2 + c\theta^3 +
#' d\theta^4 + e\theta^5 + f\theta^6
#' }
#'
#' where \eqn{\theta} is the zenith angle in radians and
#' \eqn{a,b,\dots,f} are polynomial coefficients. Polynomial degrees up to 6
#' are supported. Missing higher-order terms are assumed to be zero.
#'
#' The resulting flat-field represents the relative directional sensitivity
#' of the imaging system under spatially uniform radiance.
#'
#' @note
#' Values greater than 1.2 are considered implausible and trigger a warning.
#'
#' @param parameters numeric vector of length between 1 and 6.
#'   Polynomial coefficients ordered as \eqn{(a, b, c, d, e, f)}.
#'
#' @inheritParams add_directions
#'
#' @return Numeric [terra::SpatRaster-class] with one layer and the same
#'   geometry as `z`. Values represent relative directional sensitivity.
#'
#' @export
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' }
ff_simple_polynomial <- function(z, parameters) {

  .assert_zenith_image(z)
  .check_vector(parameters, "numeric", sign = "any")

  if (length(parameters) == 0 || length(parameters) > 6)
    stop("`parameters` must have length between 1 and 6.", call. = FALSE)

  # Pad parameters up to degree 6
  x <- rep(0, 6)
  x[seq_along(parameters)] <- parameters
  x[is.na(x)] <- 0

  theta <- z * pi / 180

  ff <- 1 +
    x[1] * theta +
    x[2] * theta^2 +
    x[3] * theta^3 +
    x[4] * theta^4 +
    x[5] * theta^5 +
    x[6] * theta^6

  # Physical validation
  vals <- terra::values(ff, mat = FALSE)

  if (any(vals[!is.na(vals)] <= 0)) {
    stop("`flat_field` must contain strictly positive values.",
         call. = FALSE)
  }

  if (any(vals[!is.na(vals)] > 1.2))
    warning("Flat-field exceeds plausible upper bound (1.2).",
         call. = FALSE)

  names(ff) <- "flat_field"
  ff
}
