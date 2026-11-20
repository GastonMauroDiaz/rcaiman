#' Validate CIE sky models
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Evaluate CIE sky coefficients by k-fold cross-validation of relative
#' radiance, keeping the supplied solar direction and zenith reference fixed.
#'
#' @details
#' Each fold is predicted by a model fitted without its observations
#' \insertCite{Kohavi1995}{rcaiman}. The solar direction and `zenith_dn` are
#' taken from `model` and kept fixed in every fit. If the original fit optimized
#' the zenith reference, its resulting value is used. This evaluates coefficient
#' fitting conditional on those two decisions, not the complete procedure used
#' to estimate the solar direction, reference, or sky samples.
#'
#' By default, rows are assigned successively to folds `1, 2, ..., k`, repeating
#' until every row has an assignment. No shuffling is performed. The original
#' row order therefore determines which observations are withheld together.
#' Supply `fold_id` to use spatial groups or a randomized partition defined in
#' the calling script. Interleaved and spatially separated folds test different
#' sampling situations; neither compensates for unobserved parts of the sky.
#'
#' The initial skies and optimization methods are supplied explicitly, rather
#' than recovered from the winning model. By default, every standard sky is
#' tried with `"Nelder-Mead"`. Use the same fitting controls as in the original
#' analysis when evaluating that procedure. The coefficients learned from the
#' full sample are not used as starting values. Custom starts should likewise
#' be chosen independently of the withheld observations.
#'
#' Each fit retains the initial-model fallback of [fit_cie_model()]. A failure
#' to obtain a usable fold model, or non-finite predictions for withheld rows,
#' stops validation and identifies that fold. With `parallel = TRUE`, fitting
#' attempts within a fold run in parallel; folds are processed sequentially.
#' Only `"SANN"` uses randomness; set a seed in the calling script if needed.
#'
#' Residuals are observed minus predicted relative radiance. For evaluation,
#' [is_outlier()] flags residuals more than `laxity` scaled MAD from their
#' median, using [stats::mad()] with its default scaling constant. This allows
#' large departures, such as cloud features to be represented by interpolation,
#' to be excluded from the regression and summary metrics. This filtering does
#' not change the samples used to fit each fold. With zero MAD or `laxity = NULL`,
#' no observations are flagged. All predictions and observations are returned
#' in their original order, including those flagged by `is_outlier`.
#'
#' The regression uses observations as the response and predictions as the
#' explanatory variable, following \insertCite{Pineiro2008;textual}{rcaiman}.
#' MSD and its components are calculated with [calc_msd()]. R-squared describes
#' association, whereas MSD and the absolute-error measures compare predictions
#' directly with observations, without correcting their slope or intercept.
#'
#' @param model list returned by [fit_cie_model()] or [optimize_sun_angles()].
#'   Its `relative_radiance` and `sun_angles` supply the observations and fixed
#'   references. Sky directions use `zenith` and `azimuth` in radians; pixel
#'   indices are not required.
#' @param k integer of length one. Number of folds, from 2 to the number of
#'   observations. Ignored when `fold_id` is supplied.
#' @param fold_id optional positive integer vector, one value per original
#'   sky-point row. Rows sharing an identifier are withheld together. At least
#'   two identifiers are required; they need not be consecutive.
#' @param laxity optional positive numeric vector of length one. Multiplier of
#'   the scaled MAD used to flag validation residuals. Default is `3`; larger
#'   values retain more observations. Use `NULL` to evaluate all observations.
#' @inheritParams fit_cie_model
#'
#' @return Plain list with:
#' \describe{
#'   \item{lm}{Observations regressed on predictions, as an [stats::lm()] object.}
#'   \item{pred, obs}{Predicted and observed relative radiance, respectively,
#'     in the original sky-point row order. Predictions are from held-out folds.}
#'   \item{r_squared}{Regression R-squared, or `NA` if predictions or observations
#'     used for evaluation are constant.}
#'   \item{rmse}{Root mean squared error.}
#'   \item{mae, medae}{Mean and median absolute error, respectively.}
#'   \item{is_outlier}{Logical vector in original row order. Flagged rows are
#'     excluded from the regression and summary metrics, but retained in
#'     `pred` and `obs`.}
#'   \item{metric}{Mean squared deviation (MSD).}
#'   \item{components}{Named vector with `SB`, `NU`, and `LC`; see [calc_msd()].}
#'   \item{fold_id}{Fold assignment of each original row.}
#'   \item{fold_models}{Fitted models, named by fold identifier and ordered by
#'     its first occurrence in `fold_id`. Each model contains its training
#'     observations, coefficients, and fitting diagnostics.}
#' }
#'
#' @seealso [fit_cie_model()], [calc_msd()]
#' @references \insertAllCited{}
#' @export
#' @examples
#' geometry <- list(model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g)
#' sun_angles <- data.frame(zenith = pi/4, azimuth = .3)
#' sky_coef <- as.numeric(cie_table[4, 1:5])
#' r <- cie_image(d, sun_angles, sky_coef) * 1000
#' names(r) <- "Blue"
#' sampling_indices <- expand.grid(col = seq(15, 65, 5), row = seq(5, 55, 5))
#' sampling_indices <- rem_unselected(sampling_indices, !is.na(r))
#' samples <- extract_dn(r, sampling_indices, use_window = FALSE)
#' samples <- add_directions(samples, d = d)
#' rr <- calc_relative_radiance(samples, "Blue", zenith_dn = 1000)
#' fit <- fit_cie_model(rr, sun_angles, std_sky_no = c(3, 4),
#'   method = "Nelder-Mead")
#' validation <- validate_cie_model(fit, k = 3, std_sky_no = c(3, 4))
#' validation$rmse
#' plot(validation$pred, validation$obs,
#'   xlab = "Predicted relative radiance", ylab = "Observed relative radiance")
#' abline(0, 1)
#'
#' \dontrun{
#' # Withhold azimuth sectors together instead of interleaving individual rows.
#' folds <- floor(fit$relative_radiance$sky_points$azimuth / (pi/2)) + 1L
#' validation <- validate_cie_model(fit, fold_id = folds,
#'   std_sky_no = c(3, 4), method = "Nelder-Mead")
#' validation$components
#' validation$fold_models[[1]]$attempts
#' }
validate_cie_model <- function(model, k = 10, fold_id = NULL,
                               custom_sky_coef = NULL, std_sky_no = NULL,
                               general_sky_type = NULL,
                               method = "Nelder-Mead", parallel = FALSE,
                               cores = NULL, logical = TRUE, leave_free = 1,
                               laxity = 3) {
  .check_vector(laxity, "numeric", 1, sign = "positive", allow_null = TRUE)
  result <- .cie_validation_predictions(model, k, fold_id,
    custom_sky_coef, std_sky_no, general_sky_type, method,
    parallel, cores, logical, leave_free)
  error <- result$obs - result$pred
  outlier <- if (is.null(laxity)) rep(FALSE, length(error)) else
    is_outlier(error, laxity = laxity)
  pred <- result$pred[!outlier]
  obs <- result$obs[!outlier]
  reg <- stats::lm(obs ~ pred)
  comparison <- calc_msd(pred, obs)
  r_squared <- if (length(unique(pred)) < 2L || length(unique(obs)) < 2L)
    NA_real_ else stats::cor(pred, obs)^2
  list(lm = reg, pred = result$pred, obs = result$obs,
    r_squared = r_squared, rmse = sqrt(unname(comparison["MSD"])),
    mae = mean(abs(obs - pred)), medae = stats::median(abs(obs - pred)),
    is_outlier = outlier, metric = unname(comparison["MSD"]),
    components = comparison[c("SB", "NU", "LC")],
    fold_id = result$fold_id, fold_models = result$fold_models)
}

# Refit each training set without using coefficients learned from held-out rows.
.cie_validation_predictions <- function(model, k, fold_id,
                                        custom_sky_coef, std_sky_no,
                                        general_sky_type, method,
                                        parallel, cores, logical, leave_free) {
  .assert_cie_model(model)
  relative_radiance <- model$relative_radiance
  if (!is.list(relative_radiance) ||
      !all(c("sky_points", "zenith_dn") %in% names(relative_radiance))) {
    stop("`model$relative_radiance` must contain `sky_points` and `zenith_dn`.",
         call. = FALSE)
  }
  samples <- relative_radiance$sky_points
  if (!is.data.frame(samples) || nrow(samples) < 2L) {
    stop("Validation requires at least two sky samples.", call. = FALSE)
  }
  .check_columns(samples, c("zenith", "azimuth", "relative_radiance"))
  for (column in c("zenith", "azimuth", "relative_radiance")) {
    .check_vector(samples[[column]], "numeric", name = column)
  }
  .check_vector(relative_radiance$zenith_dn, "numeric", 1, sign = "positive")
  .check_cie_sun_angles(model$sun_angles)
  if (any(samples$zenith < 0 | samples$zenith > pi/2) ||
      any(samples$azimuth < 0 | samples$azimuth >= 2*pi)) {
    stop("Sky directions must be in radians: zenith in [0, pi/2], azimuth in [0, 2*pi).",
         call. = FALSE)
  }
  n <- nrow(samples)
  if (is.null(fold_id)) {
    .check_vector(k, "integerish", 1, sign = "positive")
    if (k < 2L || k > n) {
      stop("`k` must be between 2 and the number of sky samples.", call. = FALSE)
    }
    fold_id <- rep(seq_len(k), length.out = n)
  } else {
    .check_vector(fold_id, "integerish", n, sign = "positive")
    if (length(unique(fold_id)) < 2L) {
      stop("`fold_id` must identify at least two folds.", call. = FALSE)
    }
  }
  groups <- unique(fold_id)
  fits <- stats::setNames(vector("list", length(groups)), as.character(groups))
  predicted <- numeric(n)
  for (j in seq_along(groups)) {
    held_out <- which(fold_id == groups[j])
    training <- relative_radiance
    training$sky_points <- samples[-held_out, , drop = FALSE]
    training$sky_points$pred <- NULL
    fits[[j]] <- tryCatch(
      fit_cie_model(training, model$sun_angles,
        custom_sky_coef = custom_sky_coef, std_sky_no = std_sky_no,
        general_sky_type = general_sky_type, optim_zenith_dn = FALSE,
        method = method, parallel = parallel, cores = cores,
        logical = logical, leave_free = leave_free),
      error = function(e) stop("Fold ", groups[j], ": ", conditionMessage(e),
                               call. = FALSE))
    predict <- .cie_cached_predictor(samples$zenith[held_out],
      samples$azimuth[held_out], model$sun_angles)
    values <- predict(fits[[j]]$coef)
    if (any(!is.finite(values))) {
      stop("Fold ", groups[j], " produced non-finite held-out predictions.",
           call. = FALSE)
    }
    predicted[held_out] <- values
  }
  list(pred = predicted, obs = samples$relative_radiance,
       fold_id = fold_id, fold_models = fits)
}
