# Design record (2026-09-04): semantic extension, not a migration of the
# complete experiment in data-raw/proyectar_pan.R. Image predictors are shared
# by calibration and prediction; the model does not construct z or a rasters.
#
# Migration record (2026-09-04, semantic consistency): interpolate_planar(),
# triangulate_planar(), fit_trend_surface(), and to_panoramic() now use x/y_down
# without snapping continuous input to cells. from_panoramic() uses the
# canonical containing-cell conversion. No interpolation engine was replaced.
#
# TODO (2026-09-04, transition): migrate higher-level callers separately.
# ff_trend_surface() still calls row_col_from_zen_azi() and constructs row/col;
# sky_via_local_models() still mixes query_points/query_indices and passes
# indices to interpolate_planar(). Do not rename indices to x/y_down as a
# compatibility workaround; use coord_from_indices() at real boundaries.
# reproject_image() already receives coord from reproject_points(), but its
# target template/mask still use the source dimensions and k/p are not passed
# to interpolation. These existing issues require their own review.
# The panoramic helpers also remain unexported, as before this migration.
#
# Correction record (2026-09-07): point transformations now remove
# relative_radius and image_angle with keep_coord() after changing the
# image frame. Recompute these predictors explicitly with add_image_polar()
# before evaluating a model in the new frame.
#
# Design record (2026-09-08, updated 2026-09-11): rcaiman defines the parameter
# contract for `model = "anisotropic_projection"`; rcaiman.registry
# operationalizes it during cross-calibration and delegates its structural
# assertion back to .assert_anisotropic_projection_parameters(). Do not
# duplicate the parameter checks there.

#' Add image-centered polar predictors
#'
#' @description
#' Enrich continuous image coordinates with relative radius and image-plane
#' angle for model fitting and prediction.
#'
#' @details
#' The image center is at `c(width / 2, height / 2)` in the `coord`
#' representation described in [y_down_from_y()]. Offsets from this center
#' retain the image axes: horizontal offsets increase rightward, and vertical
#' offsets increase downward. They are passed to [geometry::cart2pol()] in
#' that order, without converting to the ordinary `xy` representation.
#'
#' Radial distances are divided by the image semidiagonal:
#' \deqn{r_{rel} = \frac{\sqrt{(x-w/2)^2 + (y_{down}-h/2)^2}}
#'   {\sqrt{(w/2)^2 + (h/2)^2}}}
#' Thus, the center has relative radius zero and the image corners have
#' relative radius one. This normalization is not relative to the horizon.
#'
#' The image angle is measured from the rightward horizontal axis and increases
#' clockwise because the vertical image axis points downward. It is an
#' image-plane predictor, not the viewing azimuth `azi`. At the exact image
#' center, the angle is zero following the `atan2(0, 0)` convention used by
#' [geometry::cart2pol()]; this value does not identify a unique direction.
#'
#' Use this same preprocessing for calibration correspondences and for new
#' positions at which the fitted model is evaluated. Observed `zen` and `azi`
#' columns, if present, are preserved for fitting. No image raster is required.
#'
#' @note
#' Dimensions must describe the same image frame as the input coordinates.
#' Points outside that frame are retained without clipping. Convert cell
#' indices explicitly with [coord_from_indices()] when sampling cell centers.
#'
#' @param sampling_coord `data.frame` with finite numeric columns `"x"` and
#'   `"y_down"`. Additional columns are preserved. The output column names
#'   `"relative_radius"` and `"image_angle"` must not already exist.
#' @param width,height integer-like numeric vectors of length one. Positive
#'   image width and height in pixels, respectively.
#'
#' @return `data.frame` with all input columns and the following added columns:
#' \describe{
#'   \item{`relative_radius`}{Distance to the image center divided by the
#'     image semidiagonal.}
#'   \item{`image_angle`}{Image-plane angle in radians, in `[-pi, pi]`.}
#' }
#' Row order, duplicates, and zero-row inputs with the required columns are
#' preserved.
#'
#' @seealso [add_directions_from_model()], [coord_from_indices()]
#'
#' @md
#' @export
#' @examples
#' sampling_coord <- data.frame(
#'   x = c(400, 700, 400),
#'   y_down = c(100, 300, 500),
#'   id = c("top", "right", "bottom")
#' )
#' add_image_polar(sampling_coord, width = 800, height = 600)
#'
#' sampling_coord <- coord_from_indices(col = c(1, 800), row = c(1, 600))
#' add_image_polar(sampling_coord, width = 800, height = 600)
add_image_polar <- function(sampling_coord, width, height) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down"),
    reserved_columns = c("relative_radius", "image_angle")
  )
  .check_vector(width, "integerish", 1, sign = "positive")
  .check_vector(height, "integerish", 1, sign = "positive")
  for (column in c("x", "y_down")) {
    .check_vector(
      sampling_coord[[column]],
      "numeric",
      name = paste0("sampling_coord$", column)
    )
  }
  if (any(!is.finite(c(sampling_coord$x, sampling_coord$y_down)))) {
    stop("`x` and `y_down` must contain finite values.", call. = FALSE)
  }

  image_offset <- cbind(
    sampling_coord$x - width / 2,
    sampling_coord$y_down - height / 2
  )
  polar <- geometry::cart2pol(image_offset)
  semidiagonal <- sqrt((width / 2)^2 + (height / 2)^2)

  sampling_coord$relative_radius <- polar[, "r"] / semidiagonal
  sampling_coord$image_angle <- polar[, "theta"]
  sampling_coord
}


#' Add directions from a fitted image model
#'
#' @description
#' Enrich image-centered polar predictors with modeled zenith and azimuth
#' directions.
#'
#' @details
#' This function evaluates one specific model, rather than dispatching among
#' arbitrary projection families. It neither fits the coefficients nor
#' constructs the angular rasters used by [add_directions()].
#'
#' Prepare both calibration and prediction coordinates with
#' [add_image_polar()]. For relative radius \eqn{r} and image angle
#' \eqn{\phi}, define \eqn{t = \arctan(r / s)}, where \eqn{s} is
#' `atan_scale`. The model is:
#' \deqn{zen = \sum_{k=1}^{p} b_k t^k}
#' \deqn{S = u_1 \sin(\phi) + u_2 \cos(\phi)}
#' \deqn{C = v_1 \sin(\phi) + v_2 \cos(\phi)}
#' \deqn{azi = \operatorname{atan2}(S, C) \bmod 2\pi}
#' The model has no zenith intercept, so the image center satisfies
#' `zen = 0`. The vectors `zen_coef`, `azi_sin_coef`, and `azi_cos_coef`
#' contain `c(b1, ..., bp)`, `c(u1, u2)`, and `c(v1, v2)`, respectively.
#' Positions, not coefficient names, determine their interpretation.
#'
#' The azimuth coefficients must be fitted to the sine and cosine of viewing
#' azimuths in the convention of [add_directions()]. The image-plane angle is
#' only a predictor; no fixed North-up rotation is added to the predictions.
#' The same positive `atan_scale` must be used during fitting and evaluation.
#'
#' @section Fitting:
#' Enrich calibration coordinates carrying observed `zen` and `azi` with
#' [add_image_polar()], then compute `atan_radius` as
#' `atan(relative_radius / atan_scale)`. For a selected `degree`, fit the
#' following formulas with [stats::lm()]:
#' \itemize{
#'   \item `zen ~ poly(atan_radius, degree, raw = TRUE) - 1`;
#'   \item `sin(azi) ~ sin(image_angle) + cos(image_angle) - 1`;
#'   \item `cos(azi) ~ sin(image_angle) + cos(image_angle) - 1`.
#' }
#' Supply their [stats::coef()] vectors in `parameters` as `zen_coef`,
#' `azi_sin_coef`, and `azi_cos_coef`, respectively. Prepare prediction
#' positions independently with [add_image_polar()] using the same image frame.
#' Retain observed angles separately when evaluating predictions at the
#' calibration positions.
#'
#' @note
#' This is a design-stage interface. Structural input checks do not establish
#' calibration accuracy, angular coverage, or invertibility. Predicted zenith
#' angles are not clipped to the hemisphere. Jointly zero sine and cosine
#' predictions have no defined azimuth; the current arithmetic follows R's
#' `atan2(0, 0)` convention and must not be interpreted as a calibrated North.
#'
#' [reproject_points()] consumes the resulting directions and returns target
#' `x`, `y_down` coordinates for [interpolate_planar()], [triangulate_planar()],
#' or [fit_trend_surface()]. Alternatively, [to_panoramic()] produces
#' continuous coordinates and a template for planar interpolation in the
#' panoramic domain. In \pkg{rcaiman.registry}, these parameters correspond to
#' `model = "anisotropic_projection"`. The zenith component remains radial;
#' anisotropy refers to the fitted mapping from image-plane angle to azimuth.
#' Point transformations remove `relative_radius` and `image_angle` because
#' they describe the source frame. Apply [add_image_polar()] again before
#' evaluating a model in a different image frame.
#'
#' @param sampling_coord `data.frame` prepared by [add_image_polar()], with
#'   finite numeric columns `"x"`, `"y_down"`, `"relative_radius"`, and
#'   `"image_angle"`. Relative radii must be nonnegative. Additional columns
#'   are preserved. Columns `"zen"` and `"azi"` must not already exist.
#' @param parameters named list defining the model. It must contain exactly:
#'   \describe{
#'     \item{`zen_coef`}{Non-empty finite numeric vector containing the zenith
#'       polynomial coefficients in increasing order of degree, without an
#'       intercept. Its length is the polynomial degree.}
#'     \item{`azi_sin_coef`}{Finite numeric vector of length two predicting
#'       `sin(azi)` from `sin(image_angle)` and `cos(image_angle)`, in that
#'       order and without an intercept.}
#'     \item{`azi_cos_coef`}{Finite numeric vector of length two predicting
#'       `cos(azi)` from `sin(image_angle)` and `cos(image_angle)`, in that
#'       order and without an intercept.}
#'     \item{`atan_scale`}{Positive finite numeric vector of length one used in
#'       `atan(relative_radius / atan_scale)`.}
#'   }
#'
#' @return `data.frame` with all input columns and the following added columns:
#' \describe{
#'   \item{`zen`}{Predicted zenith angle in radians.}
#'   \item{`azi`}{Predicted azimuth angle in radians, in `[0, 2 * pi)`.}
#' }
#' Original coordinates, row order, duplicates, and zero-row inputs with the
#' required columns are preserved.
#'
#' @seealso [add_image_polar()], [add_directions()], [reproject_points()]
#'
#' @md
#' @export
#' @examples
#' sampling_coord <- data.frame(
#'   x = c(400, 700, 400),
#'   y_down = c(100, 300, 500),
#'   dn = c(0.2, 0.4, 0.6)
#' )
#' sampling_coord <- add_image_polar(sampling_coord, 800, 600)
#'
#' # Illustrative coefficients, not a camera calibration.
#' parameters <- list(
#'   zen_coef = 1,
#'   azi_sin_coef = c(0, -1),
#'   azi_cos_coef = c(-1, 0),
#'   atan_scale = 0.6
#' )
#' directed_coord <- add_directions_from_model(sampling_coord, parameters)
#' reproject_points(directed_coord, equidistant_coef, diameter = 100)
add_directions_from_model <- function(sampling_coord, parameters) {
  predictor_columns <- c("x", "y_down", "relative_radius", "image_angle")
  .check_columns(
    sampling_coord,
    required_columns = predictor_columns,
    reserved_columns = c("zen", "azi")
  )
  for (column in predictor_columns) {
    .check_vector(
      sampling_coord[[column]],
      "numeric",
      name = paste0("sampling_coord$", column)
    )
    if (any(!is.finite(sampling_coord[[column]]))) {
      stop("`", column, "` must contain finite values.", call. = FALSE)
    }
  }
  .check_vector(
    sampling_coord$relative_radius,
    "numeric",
    sign = "nonnegative",
    name = "sampling_coord$relative_radius"
  )
  .assert_anisotropic_projection_parameters(parameters)

  zen_coef <- parameters$zen_coef
  azi_sin_coef <- parameters$azi_sin_coef
  azi_cos_coef <- parameters$azi_cos_coef
  atan_scale <- parameters$atan_scale

  atan_radius <- atan(sampling_coord$relative_radius / atan_scale)
  sampling_coord$zen <- numeric(length(atan_radius))
  for (degree in seq_along(zen_coef)) {
    sampling_coord$zen <- sampling_coord$zen +
      zen_coef[degree] * atan_radius^degree
  }

  sin_image <- sin(sampling_coord$image_angle)
  cos_image <- cos(sampling_coord$image_angle)
  sin_azi <- azi_sin_coef[1] * sin_image + azi_sin_coef[2] * cos_image
  cos_azi <- azi_cos_coef[1] * sin_image + azi_cos_coef[2] * cos_image

  # Convention declared on 2026-09-04: coefficients fit sin(azi) and cos(azi),
  # so reconstruction uses atan2(sin_azi, cos_azi). Unlike the old experiment,
  # fitting and prediction share the same (x, y_down) preprocessing. This is
  # not a compatibility patch for its swapped axes or reversed atan2 arguments.
  # TODO (2026-09-04): validate this convention on the archived paired images;
  # determine the policy for an undefined azimuth when both predictions are 0.
  sampling_coord$azi <- atan2(sin_azi, cos_azi) %% (2 * pi)
  # Floating-point wrapping of a tiny negative angle can round up to 2 * pi.
  # Represent that endpoint as zero to honor the half-open azimuth interval.
  sampling_coord$azi[sampling_coord$azi == 2 * pi] <- 0
  sampling_coord
}
