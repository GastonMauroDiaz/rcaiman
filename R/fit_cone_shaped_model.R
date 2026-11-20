#' Fit a cone-shaped sky model
#'
#' @description
#' Fit sky digital numbers as a function of external viewing directions and
#' return parameters for [cone_shaped_image()].
#'
#' @details
#' Both models are fitted by ordinary least squares with [stats::lm()].
#' For zenith angle \eqn{\theta} and azimuth \eqn{\phi}, in radians:
#' \describe{
#'   \item{`"zenith_only"`}{\deqn{sDN = a + b\theta + c\theta^2}}
#'   \item{`"zenith_n_azimuth"`}{\deqn{sDN = a + b\theta + c\theta^2 + d\sin(\phi) + e\cos(\phi)}}
#' }
#' These are external angles, relative to the upward vertical and geographic
#' north. Supply the acquisition orientation to [direction_image()] and use
#' [add_directions()] to associate sampled values with these directions.
#' The model is independent of the image projection and optical-axis direction.
#'
#' Select suitable sky observations before fitting. Required columns must have
#' finite values; rows are not removed automatically. At least 20 observations
#' are required, and their directions must allow all coefficients to be
#' estimated. An insufficient or rank-deficient sample returns `NULL` with a
#' warning. No normalization or constraint on predicted values is imposed.
#'
#' In the model with azimuth, the sinusoidal terms do not vanish at zenith.
#' Consequently, the prediction there depends on the supplied azimuth convention.
#' A missing azimuth cannot be inferred from a zero zenith angle. The zenith-only
#' model does not use azimuth.
#'
#' See \insertCite{Diaz2018;textual}{rcaiman} for the approach underlying these
#' models. [cone_shaped_image()] evaluates the fitted parameters on an
#' image grid; [stats::predict()] applied to the returned `lm` object evaluates
#' them at directions supplied as a table.
#'
#' @param sampling_directions data frame with a numeric column `zenith` in
#'   radians and the values selected by `col_id`. For `"zenith_n_azimuth"`,
#'   a numeric column `azimuth` in radians is also required. Zenith angles must
#'   lie in \eqn{[0, \pi]} and azimuths in \eqn{[0, 2\pi)}. Additional columns,
#'   including cell indices, are ignored.
#' @param col_id character or integer vector of length one. Name or position
#'   of the column containing sampled values in `sampling_directions`.
#' @param method character vector of length one. Either `"zenith_n_azimuth"`
#'   (default), adding sinusoidal azimuth terms to a quadratic zenith model,
#'   or `"zenith_only"`, using only the quadratic zenith model.
#' @param return character vector of length one. Either `"parameters"`
#'   (default) or `"all"`.
#'
#' @return With `return = "parameters"`, a plain named list containing `method`
#'   and `coefficients`, the named coefficient vector from [stats::coef()].
#'   With `return = "all"`, a list containing `parameters` and `model`, the
#'   complete `lm` object. Use [stats::fitted()] and [stats::residuals()] for
#'   values in input row order, and `predict(model, newdata = query_directions)`
#'   for predictions from `zenith` and, when used, `azimuth` in radians.
#'   Returns `NULL` with a warning for fewer than 20 observations or a
#'   rank-deficient fit.
#'
#' @seealso [cone_shaped_image()], [extract_dn()], [add_directions()]
#' @export
#' @references \insertAllCited{}
#' @examples
#' # Illustrative sky image. With measured data, use a sky mask to select samples.
#' geometry <- list(model = "radial_projection", parameters = equidistant_coef,
#'   center_of_vision = c(40, 30), normalization_radius = 30,
#'   normalization_angle = 90)
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 60)
#' r <- 1000 + 200 * (d$zenith * pi / 180)^2 +
#'   50 * sin(d$azimuth * pi / 180) + 10 * cos(2 * d$azimuth * pi / 180)
#' names(r) <- "Blue"
#' m <- !is.na(d$zenith) & !is.na(d$azimuth) & d$zenith < 90
#' sampling_indices <- expand.grid(col = seq(5, 75, 5), row = seq(5, 55, 5))
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' samples <- extract_dn(r, sampling_indices, use_window = FALSE)
#' samples <- add_directions(samples, d = d)
#' fit <- fit_cone_shaped_model(samples, col_id = "Blue", return = "all")
#' summary(fit$model)
#' sky <- cone_shaped_image(d, fit$parameters, m = m)
#' plot(sky)
#' query_directions <- sky_grid_centers(15)
#' head(predict(fit$model, newdata = query_directions))
fit_cone_shaped_model <- function(sampling_directions,
                                  col_id,
                                  method = "zenith_n_azimuth",
                                  return = "parameters") {
  .assert_choice(method, c("zenith_only", "zenith_n_azimuth"))
  .assert_choice(return, c("parameters", "all"))
  columns <- if (method == "zenith_only") "zenith" else c("zenith", "azimuth")
  .check_columns(sampling_directions, columns)
  if (is.numeric(col_id)) .check_vector(col_id, "integerish", 1, sign = "positive")
  col_id <- .resolve_required_column(col_id, sampling_directions)
  data <- sampling_directions[columns]
  for (column in columns) {
    angle <- data[[column]]
    .check_vector(angle, "numeric", name = column)
    invalid <- if (column == "zenith") angle < 0 | angle > pi else
      angle < 0 | angle >= 2 * pi
    if (any(invalid)) {
      stop("`", column, "` must be in ",
           if (column == "zenith") "[0, pi]" else "[0, 2*pi)",
           " radians.", call. = FALSE)
    }
  }
  data$observed <- sampling_directions[[col_id]]
  .check_vector(data$observed, "numeric", name = col_id)
  if (nrow(data) < 20) {
    warning("At least 20 observations are required to fit the model.", call. = FALSE)
    return(NULL)
  }
  formula <- if (method == "zenith_only") {
    observed ~ zenith + I(zenith^2)
  } else {
    observed ~ zenith + I(zenith^2) + sin(azimuth) + cos(azimuth)
  }
  model <- stats::lm(formula, data = data, na.action = stats::na.fail)
  coefficients <- stats::coef(model)
  if (model$rank < length(coefficients) || any(!is.finite(coefficients))) {
    warning("Sample directions do not support a full-rank model.", call. = FALSE)
    return(NULL)
  }
  parameters <- list(method = method, coefficients = coefficients)
  if (return == "parameters") return(parameters)
  list(parameters = parameters, model = model)
}
