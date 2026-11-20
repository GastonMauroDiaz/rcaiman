#' Fit a polynomial trend surface
#'
#' @description
#' Fit a polynomial surface to observations in image coordinates and return
#' compact parameters for evaluating it on a raster grid.
#'
#' @details
#' The surface is fitted with [spatial::surf.ls()] using `x` and `y_down`
#' directly, following [coord_from_indices()]. Values may be positive, zero,
#' or negative. No normalization is performed and the constant term is
#' fitted freely.
#'
#' A degree `np` surface has `(np + 1) * (np + 2) / 2` coefficients. Sample
#' positions must provide at least that many independent constraints. The
#' coordinate ranges retained in the parameters scale the polynomial basis;
#' they are not image dimensions and must not be replaced during evaluation.
#'
#' [trend_surface_image()] evaluates the parameters on a grid. To restrict
#' predictions to the convex hull of the samples, apply [convex_hull_mask()]
#' separately. The parameters contain neither observations nor their hull.
#'
#' For an application of this method to radiometric preprocessing, see the
#' "Radiometric preprocessing" section in [trend_surface_image()].
#'
#' @param sampling_coord `data.frame` with numeric columns `x`, `y_down`, and
#'   the values selected by `col_id`. Coordinates and values must be finite.
#'   Remove missing observations before fitting.
#' @param col_id character or numeric vector of length one. Name or position
#'   of the column to model in `sampling_coord`.
#' @param np integer vector of length one. Polynomial degree, from 1 to 6.
#' @param return character vector of length one. Either `"parameters"`
#'   (default) or `"all"`.
#'
#' @return With `return = "parameters"`, a plain named list containing:
#'   \describe{
#'     \item{`np`}{polynomial degree.}
#'     \item{`beta`}{fitted coefficients in the basis used by `spatial`.}
#'     \item{`rx`, `ry`}{ranges of `x` and `y_down` used during fitting.}
#'   }
#'   With `return = "all"`, a list containing `parameters`, `model` (the
#'   complete `trls` object returned by [spatial::surf.ls()]), and `r2`.
#'   Use [stats::fitted()] and [stats::residuals()] on `model` for fitted values
#'   and residuals in input row order. For predictions at other positions,
#'   use `predict(model, x, y_down)`; the model's `y` coordinate is `y_down`.
#'   `r2` is `NA` for constant observations.
#'
#' @seealso [trend_surface_image()], [convex_hull_mask()]
#' @export
#' @examples
#' samples <- expand.grid(x = c(2, 6, 10, 14, 18), y_down = c(2, 5, 8, 10))
#' samples$dn <- 10 + 2 * samples$x - 3 * samples$y_down +
#'   0.1 * samples$x * samples$y_down
#' fit <- fit_trend_surface(samples, np = 2, return = "all")
#' range(residuals(fit$model))
#' r <- set_rcaiman_geometry(terra::rast(nrows = 12, ncols = 20))
#' m <- convex_hull_mask(samples, r)
#' surface <- trend_surface_image(r, fit$parameters, m = m)
#' plot(surface)
fit_trend_surface <- function(sampling_coord,
                                 col_id = "dn",
                                 np = 6,
                                 return = "parameters") {
  .check_columns(sampling_coord, required_columns = c("x", "y_down"))
  .check_vector(np, "integerish", 1, sign = "positive")
  if (np > 6) stop("`np` must lie within [1, 6].", call. = FALSE)
  .assert_choice(return, c("parameters", "all"))
  if (is.numeric(col_id)) {
    .check_vector(col_id, "integerish", 1, sign = "positive")
    if (col_id > ncol(sampling_coord)) {
      stop("`col_id` exceeds the number of columns.", call. = FALSE)
    }
    col_id <- names(sampling_coord)[col_id]
  } else {
    .check_vector(col_id, "character", 1)
  }
  .check_columns(sampling_coord, required_columns = col_id)
  x <- sampling_coord$x
  y_down <- sampling_coord$y_down
  observed <- sampling_coord[[col_id]]
  .check_vector(x, "numeric")
  .check_vector(y_down, "numeric")
  .check_vector(observed, "numeric")
  npar <- (np + 1) * (np + 2) / 2
  if (length(observed) < npar) {
    stop("Degree ", np, " requires at least ", npar, " observations.",
         call. = FALSE)
  }
  if (diff(range(x)) == 0 || diff(range(y_down)) == 0) {
    stop("Both coordinate ranges must have positive width.", call. = FALSE)
  }
  model <- spatial::surf.ls(np, x, y_down, observed)
  # surf.ls() does not propagate its native solver's rank-failure flag.
  if (qr(matrix(model$f, nrow = length(observed)))$rank < npar ||
      any(!is.finite(model$beta))) {
    stop("Sample positions do not support a full-rank surface of degree ",
         np, ".", call. = FALSE)
  }
  parameters <- model[c("np", "beta", "rx", "ry")]
  if (return == "parameters") return(parameters)

  residuals <- stats::residuals(model)
  ss_tot <- sum((observed - mean(observed))^2)
  r2 <- if (ss_tot == 0) NA_real_ else 1 - sum(residuals^2) / ss_tot
  list(parameters = parameters, model = model, r2 = r2)
}

#' @rdname internal-validation
#' @export
.check_trend_surface_parameters <- function(parameters) {
  required <- c("np", "beta", "rx", "ry")
  if (!is.list(parameters) || !all(required %in% names(parameters)) ||
      anyDuplicated(names(parameters))) {
    stop("`parameters` must be a named list with `np`, `beta`, `rx`, and `ry`.",
         call. = FALSE)
  }
  .check_vector(parameters$np, "integerish", 1, sign = "positive")
  if (parameters$np > 6) stop("`np` must lie within [1, 6].", call. = FALSE)
  npar <- (parameters$np + 1) * (parameters$np + 2) / 2
  .check_vector(parameters$beta, "numeric", npar)
  .check_vector(parameters$rx, "numeric", 2)
  .check_vector(parameters$ry, "numeric", 2)
  if (diff(parameters$rx) <= 0 || diff(parameters$ry) <= 0) {
    stop("`rx` and `ry` must be increasing ranges.", call. = FALSE)
  }
  invisible(TRUE)
}

.predict_trend_surface <- function(parameters, x, y_down) {
  # Rebuild only the state used by spatial's prediction method.
  model <- structure(parameters[c("np", "beta", "rx", "ry")], class = "trls")
  # stats::predict(model, x, y_down)
  spatial::predict.trls(model, x, y_down)
}
