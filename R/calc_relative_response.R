#' Evaluate relative response and construct a flat-field image
#'
#' @description
#' Evaluate an angular vignetting model at off-axis angles with
#' `calc_relative_response()`, or construct its flat-field image with
#' `flat_field_image()`.
#'
#' @details
#' ### Relative response and correction
#'
#' Relative response describes the directional sensitivity of the imaging
#' system under spatially uniform radiance. Calibration observations are
#' normalized to the response at the center of vision, on the optical axis.
#' This reference is independent of acquisition orientation and differs from
#' the zenith radiance used by [calc_relative_radiance()] for CIE modeling.
#'
#' Evaluating the model on `g[["off_axis"]]` from [set_internal_geometry()]
#' produces a flat-field image. Divide the measured image by that flat-field
#' with [apply_flat_field()] to correct its relative response, before
#' reprojection or other interpolation of the measurements.
#'
#' ### Simple polynomial
#'
#' With `model = "simple_polynomial"`, the constant term is fixed at one:
#'
#' \deqn{f_v(\theta') = 1 + \sum_{j=1}^{p} b_j\theta'^j.}
#'
#' Here \eqn{\theta'} is the off-axis angle in radians. `parameters` contains
#' \eqn{b_1, \ldots, b_p}, in increasing order of degree, for degrees 1 to 6.
#' The model preserves the optical-axis reference exactly: \eqn{f_v(0) = 1}.
#'
#' ### Free-form polynomial
#'
#' With `model = "free_form_polynomial"`, the constant term is fitted:
#'
#' \deqn{f_v(\theta') = \sum_{j=0}^{p} b_j\theta'^j.}
#'
#' `parameters` contains \eqn{b_0, \ldots, b_p}, for degrees 0 to 6, and
#' \eqn{\theta'} is in radians. Calibration observations use the same
#' optical-axis normalization as for the simple polynomial, but
#' \eqn{f_v(0) = b_0} need not equal one. Normalizing observations and
#' constraining a fitted model are separate choices. Predictions are not
#' renormalized. For either polynomial, supply only the coefficients through
#' the fitted degree; absent higher-order terms are zero.
#'
#' ### Lang (2010)
#'
#' With `model = "lang2010"`, Eq. (5) from
#' \insertCite{Lang2010;textual}{rcaiman} is evaluated for the blue band
#' reported in that publication:
#'
#' \deqn{f_v(\theta') = 0.622 - \frac{1.504\theta'}{A^{4.553}} +
#' 0.0539\ln((22.87-A)(91-\theta')).}
#'
#' In this expression \eqn{\theta'} is in degrees, between 0 and 90, and
#' \eqn{A} is the acquisition aperture f-number. Unit conversion is handled
#' internally. The empirical coefficients are built in: leave `parameters`
#' as `NULL` and supply `aperture`. The expression refers to the center of
#' vision but does not equal exactly one there; it is not renormalized.
#'
#' ### Inspecting and applying a model
#'
#' `calc_relative_response()` returns predictions, including zero or negative
#' values, so that the fitted curve can be inspected. Use angles within the
#' range supported by the calibration. `flat_field_image()` requires finite,
#' strictly positive predictions at non-missing angles. For the simple
#' polynomial, values above 1.2 trigger a plausibility warning.
#'
#' For a spatial response surface rather than an angular model, see the
#' "Radiometric preprocessing" section in [trend_surface_image()].
#'
#' @param off_axis numeric vector of off-axis angles for
#'   `calc_relative_response()`, in the units declared by `units`.
#'   For `flat_field_image()`, a single-layer [terra::SpatRaster-class] of
#'   off-axis angles in degrees, normally `g[["off_axis"]]`.
#'   Angles must lie within `[0, 180]` degrees (`[0, pi]` radians), or within
#'   `[0, 90]` degrees for `lang2010`. Missing values are retained.
#' @param model character vector of length one: `"simple_polynomial"`,
#'   `"free_form_polynomial"`, or `"lang2010"`.
#' @param parameters finite numeric vector of polynomial coefficients:
#'   length 1 to 6 without an intercept for `simple_polynomial`, or length
#'   1 to 7 including the intercept for `free_form_polynomial`.
#'   Must be `NULL` for `lang2010`.
#' @param units character vector of length one, `"degree"` or `"radian"`.
#'   Declares the input vector's units. Polynomial coefficients always apply
#'   to radians, irrespective of this argument. Tables from [add_directions()]
#'   contain radians; declare `units = "radian"` when evaluating them.
#' @param aperture numeric vector of length one. Acquisition aperture f-number
#'   for `lang2010`, greater than zero and less than 22.87. Leave as `NULL`
#'   for polynomial models, whose coefficients already describe the aperture
#'   at which they were calibrated.
#'
#' @return `calc_relative_response()` returns a numeric vector of dimensionless
#'   responses with the input's length, order, and names.
#'   `flat_field_image()` returns a single-layer [terra::SpatRaster-class]
#'   named `flat_field`, preserving the input grid and missing values.
#'
#' @references \insertAllCited{}
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(40, 30),
#'   parameters = 2 / pi,
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' # Illustrative coefficients; use those fitted to your calibration.
#' vignetting_coef <- c(-0.1, -0.05)
#' points <- add_directions(data.frame(col = c(40, 45), row = c(30, 35)), g)
#' calc_relative_response(points$off_axis, "simple_polynomial",
#'                        vignetting_coef, units = "radian")
#' angle <- seq(0, 90, length.out = 100)
#' response <- calc_relative_response(angle, "simple_polynomial", vignetting_coef)
#' plot(angle, response, type = "l", xlab = "Off-axis angle (deg)",
#'      ylab = "Relative response")
#' ff <- flat_field_image(g[["off_axis"]], "simple_polynomial", vignetting_coef)
#' plot(ff)
#' measured <- ff * 100
#' corrected <- apply_flat_field(measured, ff)
#' calc_relative_response(c(0, 30, 60, 90), "lang2010", aperture = 5)
calc_relative_response <- function(off_axis, model, parameters = NULL,
                                   units = "degree", aperture = NULL) {
  .check_relative_response_model(model, parameters, aperture)
  .check_vector(off_axis, "numeric", allow_na = TRUE, name = "off_axis")
  .check_vector(units, "character", 1)
  .assert_choice(units, c("degree", "radian"))
  limit <- if (model == "lang2010") 90 else 180
  if (units == "radian") limit <- limit * pi / 180
  if (any(off_axis < 0 | off_axis > limit, na.rm = TRUE))
    stop("`off_axis` must lie within [0, ", limit, "] in ", units,
         " units for this model.", call. = FALSE)

  if (model == "lang2010") {
    angle <- if (units == "radian") off_axis * 180 / pi else off_axis
    response <- 0.622 - 1.504 * angle / aperture^4.553 +
      0.0539 * log((22.87 - aperture) * (91 - angle))
  } else {
    angle <- if (units == "degree") off_axis * pi / 180 else off_axis
    coef <- if (model == "simple_polynomial") c(1, parameters) else parameters
    response <- rep(0, length(angle))
    for (i in rev(seq_along(coef))) response <- response * angle + coef[i]
  }
  names(response) <- names(off_axis)
  response
}

#' @rdname calc_relative_response
#' @export
flat_field_image <- function(off_axis, model, parameters = NULL, aperture = NULL) {
  .assert_single_layer(off_axis)
  angles <- terra::values(off_axis, mat = FALSE)
  response <- calc_relative_response(angles,
                                    model, parameters, aperture = aperture)
  valid <- !is.na(angles)
  if (any(!is.finite(response[valid]) | response[valid] <= 0))
    stop("`flat_field` must contain finite, strictly positive values.", call. = FALSE)
  if (model == "simple_polynomial" && any(response[valid] > 1.2))
    warning("Flat-field exceeds plausible upper bound (1.2).", call. = FALSE)
  ff <- terra::setValues(off_axis, response)
  names(ff) <- "flat_field"
  ff
}

.check_relative_response_model <- function(model, parameters, aperture) {
  .check_vector(model, "character", 1)
  .assert_choice(model, c("simple_polynomial", "free_form_polynomial", "lang2010"))
  if (model == "lang2010") {
    if (!is.null(parameters))
      stop("`parameters` must be NULL for lang2010.", call. = FALSE)
    .check_vector(aperture, "numeric", 1, sign = "positive")
    if (aperture >= 22.87)
      stop("`aperture` must be less than 22.87.", call. = FALSE)
  } else {
    if (!is.null(aperture))
      stop("`aperture` must be NULL for polynomial models.", call. = FALSE)
    .check_vector(parameters, "numeric")
    limit <- if (model == "simple_polynomial") 6 else 7
    if (length(parameters) == 0 || length(parameters) > limit)
      stop("`parameters` must have length between 1 and ", limit, ".", call. = FALSE)
  }
  invisible(TRUE)
}
