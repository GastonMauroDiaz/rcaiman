#' Fit a CIE sky model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Fit the CIE General Sky model to relative radiance sampled at external sky
#' directions. Compare optimized models with the selected initial models and
#' return the candidate with the smallest mean squared deviation.
#'
#' @details
#' The approach is based on \insertCite{Lang2010;textual}{rcaiman}. Sampled DN
#' should be linearly related to incident radiance and corrected for flat-field
#' response before normalization with [calc_relative_radiance()]. As a compromise
#' for gamma-encoded images, see [invert_gamma_correction()].
#'
#' Every initial coefficient set is tried with every requested optimization
#' method. Initial models also participate without optimization. They remain useful candidates when optimization fails or
#' does not improve the fit. Ties favor the first candidate evaluated.
#' `"SANN"` is stochastic; use `set.seed()` outside this function for reproducibility.
#'
#' With `parallel = TRUE`, independent start/method combinations are distributed
#' among worker processes. Results and ties retain their declared order. Worker
#' creation has a cost, so short fits may be faster sequentially.
#'
#' Optimization minimizes RMSE divided by the mean input relative radiance,
#' plus penalties. The existing penalties discourage positive `b`, negative `e`,
#' negative predictions on a hemispherical Fibonacci distribution, and sampled
#' predictions brighter than the prediction towards the sun (the latter penalty
#' is scaled by `abs(c)`).
#' When the sun is below the horizon, use the horizon at the solar azimuth for
#' this comparison. This extends the heuristic without evaluating sky radiance
#' below the hemisphere; the actual solar direction is still used in the CIE
#' indicatrix. For a sun above the horizon, the comparison is unchanged.
#' These are numerical penalties, not hard constraints or a guarantee of physical
#' validity between evaluated directions.
#'
#' The negative-prediction check uses 826 directions with approximately 5 deg
#' spacing. This density balances coverage against computation time. The penalty
#' is the squared sum of negative predictions scaled by `1e10`.
#'
#' Candidates are compared by unpenalized mean squared deviation (MSD).
#' [calc_msd()] computes this metric and its decomposition, returned in
#' `components`. See that function for the equations and their interpretation.
#' This comparison describes the fitted sample, not independent validation.
#'
#' With `optim_zenith_dn = TRUE`, the sampled DN are recovered as
#' `relative_radiance * zenith_dn` and held fixed while five sky coefficients
#' and the zenith reference are fitted jointly. The reference is represented in
#' units of the sample standard deviation of these DN. Let q0 be the initial
#' standardized reference, q the trial value, and t the smallest sampled zenith
#' angle divided by `pi/2`. The additional penalty is `max(0, abs(q - q0) - t)`.
#' Thus, more distant samples allow a wider unpenalized change. For constant
#' sampled DN, the reference is kept fixed with a warning. Candidate MSD uses
#' the relative observations normalized by that candidate's reference, so
#' changing the reference also changes the comparison scale.
#'
#' Errors from individual optimization attempts are recorded in `attempts`.
#' Finite results with nonzero convergence codes remain candidates and are
#' identified there. The function warns if all attempts fail to yield usable
#' results or if the selected optimized result did not converge. Initial models
#' are identified by `source`, with `opt_result = NULL` and `method = NA`.
#'
#' @param relative_radiance list, typically from [calc_relative_radiance()], with
#'   positive finite `zenith_dn` and a nonempty data frame `sky_points` containing
#'   finite numeric columns `zenith`, `azimuth`, and `relative_radiance`. Directions are external,
#'   in radians: zenith in `[0, pi/2]` and azimuth in `[0, 2*pi)`. Relative radiance
#'   must have a positive mean. Additional columns are preserved; pixel indices
#'   are not required. Select suitable sky samples before fitting; no rows are
#'   removed automatically.
#' @param sun_angles data frame of one row with finite numeric `zenith` and
#'   `azimuth` columns in radians, using the same external reference as the sky
#'   samples. Follows the contract of [cie_image()]: zenith is nonnegative
#'   (values above `96 * pi/180` produce a warning) and azimuth is in `[0, 2*pi)`.
#'   Additional columns are preserved and ignored during fitting.
#' @param custom_sky_coef finite numeric vector of length five or nonempty numeric
#'   matrix with five columns, in the order `a`, `b`, `c`, `d`, `e`. Custom starting
#'   coefficients, used instead of standard skies. Do not combine with standard
#'   sky selectors.
#' @param std_sky_no integer vector with standard sky numbers from 1 to 15; see
#'   [cie_table]. Takes precedence over `general_sky_type`. With neither selector,
#'   all standard skies are used.
#' @param general_sky_type character vector of length one or two, containing
#'   `"Overcast"`, `"Partly cloudy"`, or `"Clear"`. Selects initial skies from
#'   either supplied type in [cie_table], for example
#'   `c("Overcast", "Partly cloudy")` when the type is uncertain. Does not
#'   constrain the optimized coefficients to those types. With `NULL` and
#'   no `std_sky_no`, all standard skies are used.
#' @param optim_zenith_dn logical of length one. Whether to optimize the supplied
#'   zenith reference along with the sky coefficients. See *Details*.
#' @param method character vector of optimization methods passed to [stats::optim()].
#'   Supports `"Nelder-Mead"`, `"BFGS"`, `"CG"`, `"SANN"`, and `"L-BFGS-B"`.
#'
#' @param parallel logical of length one. If `TRUE`, distribute optimization
#'   attempts among worker processes. With `FALSE` (default), run attempts
#'   sequentially. SANN uses an isolated process even in sequential mode; see
#'   *Details*. The number of workers is capped at the number of attempts.
#' @inheritParams apply_locally
#'
#' @return Plain list with:
#' \describe{
#'   \item{`relative_radiance`}{Input relative radiance with the selected `zenith_dn`, observations
#'     renormalized if needed, and a `pred` column added to `sky_points`. Original
#'     columns, including any sampled DN, are preserved.}
#'   \item{`coef`}{Named vector of five coefficients for [cie_image()].}
#'   \item{`sun_angles`}{The input solar-direction table, in radians.}
#'   \item{`opt_result`}{Selected [stats::optim()] result, or `NULL` for an initial model.}
#'   \item{`method`}{Selected optimization method, or `NA` for an initial model.}
#'   \item{`start`}{Starting coefficients, followed by the standardized initial
#'     zenith reference when it is optimized.}
#'   \item{`metric`}{MSD between predicted and relative observed values.}
#'   \item{`components`}{Named numeric vector with `SB`, `NU`, and `LC`.}
#'   \item{`source`}{`"standard"`, `"custom"`, or `"optimized"`.}
#'   \item{`attempts`}{Table with initial-set index `start_id`, `source`, `method`,
#'     `convergence`, `metric`, `objective` (including penalties), `message`, and
#'     `selected`. Failed or numerically unusable candidates have `metric = Inf`.}
#' }
#'
#' @section Acquisition and sampling:
#' Use [direction_image()] to describe acquisition orientation and [add_directions()]
#' to associate sampled DN with external angles. The sun must use the same
#' external reference. The image need not include zenith: its radiometric
#' reference can be supplied separately to [calc_relative_radiance()].
#' Point files from rcaiman-projects contain continuous `x` and `y_down` coordinates;
#' use [add_indices()] before extracting DN and raster directions.
#' See *Details* for the fitting procedure and optional reference adjustment.
#'
#' @note `"SANN"` in `method` uses randomness; set a seed in the calling script for
#' repeatability. This function does not set a seed or change the caller's
#' random-number generator configuration.
#'
#' @seealso [cie_image()], [calc_relative_radiance()], [cie_table]
#' @references \insertAllCited{}
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 60)
#' sun_angles <- data.frame(zenith = 45 * pi/180, azimuth = 10 * pi/180)
#' sky_coef <- as.numeric(cie_table[4, 1:5])
#' r <- cie_image(d, sun_angles, sky_coef) * 1000
#' names(r) <- "Blue"
#' m <- !is.na(r)
#' sampling_indices <- expand.grid(col = seq(5, 75, 5), row = seq(5, 55, 5))
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' samples <- extract_dn(r, sampling_indices, use_window = FALSE)
#' samples <- add_directions(samples, d = d)
#' rr <- calc_relative_radiance(samples, "Blue", zenith_dn = 1000)
#' model <- fit_cie_model(rr, sun_angles, std_sky_no = c(3, 4), method = "Nelder-Mead")
#' model$components
#' model$attempts
#' sky <- cie_image(d, model$sun_angles, model$coef) * model$relative_radiance$zenith_dn
#' plot(sky)
#' plot(model$relative_radiance$sky_points$relative_radiance, model$relative_radiance$sky_points$pred)
#' abline(0, 1)
#'
#' \dontrun{
#' # Continue the measured-image preprocessing example in calc_relative_radiance().
#' # r is flat-field corrected; d describes the acquisition orientation.
#' path <- system.file("external/sun_disk_mark.csv",
#'                     package = "rcaiman")
#' sun_indices <- read_points(path) %>% add_indices()
#' sun_angles <- add_directions(sun_indices, d = d)
#'
#' set.seed(7)
#' model <- fit_cie_model(
#'   rr,
#'   sun_angles,
#'   general_sky_type = c("Overcast" , "Partly cloudy"),
#'   parallel = TRUE
#' )
#' sky <- cie_image(d, model$sun_angles, model$coef) * model$relative_radiance$zenith_dn
#' plot(sky)
#' plot(r / sky)
#' }
fit_cie_model <- function(relative_radiance, sun_angles,
                          custom_sky_coef = NULL,
                          std_sky_no = NULL,
                          general_sky_type = NULL,
                          optim_zenith_dn = FALSE,
                          method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
                          parallel = FALSE, cores = NULL,
                          logical = TRUE, leave_free = 1) {
  if (!is.list(relative_radiance) ||
      !all(c("sky_points", "zenith_dn") %in% names(relative_radiance))) {
    stop("`relative_radiance` must contain `sky_points` and `zenith_dn`.", call. = FALSE)
  }
  rr <- relative_radiance
  .check_vector(rr$zenith_dn, "numeric", 1, sign = "positive")
  if (!is.data.frame(rr$sky_points) || !nrow(rr$sky_points)) {
    stop("`sky_points` must be a nonempty data frame.", call. = FALSE)
  }
  .check_columns(rr$sky_points, c("zenith", "azimuth", "relative_radiance"))
  for (column in c("zenith", "azimuth", "relative_radiance")) {
    .check_vector(rr$sky_points[[column]], "numeric", name = column)
  }
  Zp <- rr$sky_points$zenith
  AzP <- rr$sky_points$azimuth
  observed <- rr$sky_points$relative_radiance
  if (any(Zp < 0 | Zp > pi/2) || any(AzP < 0 | AzP >= 2*pi)) {
    stop("Sky directions must be in radians: zenith in [0, pi/2], azimuth in [0, 2*pi).", call. = FALSE)
  }
  if (!is.finite(mean(observed)) || mean(observed) <= 0) {
    stop("Mean relative radiance must be finite and positive.", call. = FALSE)
  }
  .check_cie_sun_angles(sun_angles)
  .check_vector(optim_zenith_dn, "logical", 1)
  .check_cie_general_sky_type(general_sky_type)
  .check_vector(method, "character")
  if (!length(method) || any(!method %in% c("Nelder-Mead", "BFGS", "CG", "SANN", "L-BFGS-B"))) {
    stop("`method` must contain supported multivariate stats::optim methods.", call. = FALSE)
  }
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")
  method <- unique(method)
  if (!is.null(custom_sky_coef)) {
    if (!is.null(std_sky_no) || !is.null(general_sky_type)) {
      stop("Do not combine `custom_sky_coef` with standard sky selectors.", call. = FALSE)
    }
    if (is.numeric(custom_sky_coef) && is.null(dim(custom_sky_coef)) && length(custom_sky_coef) == 5L) {
      custom_sky_coef <- matrix(custom_sky_coef, nrow = 1)
    }
    if (!is.matrix(custom_sky_coef) || !is.numeric(custom_sky_coef) ||
        ncol(custom_sky_coef) != 5L || !nrow(custom_sky_coef) || any(!is.finite(custom_sky_coef))) {
      stop("`custom_sky_coef` must be a finite numeric vector of length five or a nonempty five-column matrix.", call. = FALSE)
    }
    starts <- custom_sky_coef
    initial_source <- "custom"
  } else {
    skies <- rcaiman::cie_table
    if (!is.null(std_sky_no)) {
      .check_vector(std_sky_no, "integerish", sign = "positive")
      if (!length(std_sky_no) || any(std_sky_no > 15)) stop("`std_sky_no` must contain numbers from 1 to 15.", call. = FALSE)
      skies <- skies[unique(std_sky_no), , drop = FALSE]
    } else if (!is.null(general_sky_type)) {
      skies <- skies[skies$general_sky_type %in% general_sky_type, , drop = FALSE]
    }
    starts <- as.matrix(skies[, 1:5])
    initial_source <- "standard"
  }

  dn <- observed * rr$zenith_dn
  if (optim_zenith_dn) {
    sigma_dn <- stats::sd(dn)
    if (!is.finite(sigma_dn) || sigma_dn == 0) {
      warning("Sampled DN have no usable standard deviation; keeping `zenith_dn` fixed.", call. = FALSE)
      optim_zenith_dn <- FALSE
    } else {
      mu_dn <- mean(dn)
      reference_score <- (rr$zenith_dn - mu_dn) / sigma_dn
      reference_tolerance <- min(Zp) / (pi/2)
    }
  }
  reference_dn <- function(params) {
    if (optim_zenith_dn) params[6] * sigma_dn + mu_dn else rr$zenith_dn
  }
  vault <- external_directions_from_internal(fibonacci_directions(5))
  predict_all <- .cie_cached_predictor(
    c(Zp, vault$zenith, min(sun_angles$zenith, pi/2)),
    c(AzP, vault$azimuth, sun_angles$azimuth), sun_angles)
  sample_index <- seq_along(observed)
  vault_index <- length(observed) + seq_len(nrow(vault))
  sun_index <- length(observed) + nrow(vault) + 1L
  predict_sky <- function(params) predict_all(params)[sample_index]
  adjust_reference <- optim_zenith_dn
  initial_reference <- rr$zenith_dn
  mean_observed <- mean(observed)
  fcost <- function(params) {
    reference <- if (adjust_reference) params[6] * sigma_dn + mu_dn else initial_reference
    if (any(!is.finite(params)) || !is.finite(reference) || reference <= 0) return(Inf)
    values <- predict_all(params)
    if (any(!is.finite(values))) return(Inf)
    pred <- values[sample_index]
    sky <- values[vault_index]
    sun <- values[sun_index]
    penalty <- 1e10 * (max(0, params[2]) + max(0, -params[5]) +
      abs(params[3]) * max(0, max(pred) - sun) + sum(pmin(sky, 0))^2)
    if (adjust_reference) {
      penalty <- penalty + max(0, abs(reference_score - params[6]) - reference_tolerance)
    }
    y <- if (adjust_reference) dn / reference else observed
    unname(sqrt(mean((pred - y)^2)) / mean_observed + penalty)
  }
  # Send only numeric fitting data to workers, not the caller's environment.
  environment(fcost) <- list2env(list(predict_all = predict_all,
    sample_index = sample_index, vault_index = vault_index, sun_index = sun_index,
    adjust_reference = adjust_reference, initial_reference = initial_reference,
    observed = observed, dn = dn, mean_observed = mean_observed,
    sigma_dn = if (optim_zenith_dn) sigma_dn else NULL,
    mu_dn = if (optim_zenith_dn) mu_dn else NULL,
    reference_score = if (optim_zenith_dn) reference_score else NULL,
    reference_tolerance = if (optim_zenith_dn) reference_tolerance else NULL),
    parent = baseenv())

  evaluate <- function(params, start, source, method = NA_character_, fit = NULL) {
    reference <- if (source == "optimized") reference_dn(params) else rr$zenith_dn
    if (any(!is.finite(params)) || !is.finite(reference) || reference <= 0) return(NULL)
    pred <- predict_sky(params)
    if (any(!is.finite(pred))) return(NULL)
    result_rr <- rr
    result_rr$zenith_dn <- unname(reference)
    if (optim_zenith_dn && source == "optimized") result_rr$sky_points$relative_radiance <- dn / reference
    result_rr$sky_points$pred <- pred
    comparison <- calc_msd(pred, result_rr$sky_points$relative_radiance)
    components <- comparison[c("SB", "NU", "LC")]
    metric <- unname(comparison["MSD"])
    if (!is.finite(metric)) return(NULL)
    list(relative_radiance = result_rr, opt_result = fit, coef = stats::setNames(params[1:5], letters[1:5]),
      sun_angles = sun_angles, method = method, start = start, metric = metric,
      components = components, source = source)
  }
  models <- list()
  attempts <- list()
  record <- function(candidate, start_id, source, method, convergence, message, objective) {
    j <- length(attempts) + 1L
    models[j] <<- list(candidate)
    attempts[[j]] <<- data.frame(start_id = start_id, source = source, method = method,
      convergence = convergence, metric = if (is.null(candidate)) Inf else candidate$metric,
      objective = objective, message = message, stringsAsFactors = FALSE)
  }
  start_vectors <- lapply(seq_len(nrow(starts)), function(i) {
    start <- stats::setNames(as.numeric(starts[i, ]), letters[1:5])
    if (optim_zenith_dn) start <- c(start, zenith_reference = reference_score)
    start
  })
  jobs <- unlist(lapply(start_vectors, function(start) {
    lapply(method, function(m) list(start = start, method = m, seed = NULL))
  }), recursive = FALSE)
  workers <- if (parallel) min(.cores(cores, logical, leave_free), length(jobs)) else 1L
  stochastic <- which(vapply(jobs, function(job) job$method == "SANN", TRUE))
  runner <- .cie_run_optim
  environment(runner) <- baseenv()
  if (workers > 1L || length(stochastic)) {
    # A single isolated worker also protects the caller's RNG in sequential SANN.
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (length(stochastic)) {
      initial_seed <- sample.int(.Machine$integer.max, 1L)
      make_streams <- function(seed, n) {
        set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
        streams <- vector("list", n)
        streams[[1]] <- get(".Random.seed", envir = .GlobalEnv)
        if (n > 1L) for (j in 2:n) streams[[j]] <- parallel::nextRNGStream(streams[[j-1L]])
        streams
      }
      environment(make_streams) <- baseenv()
      seed_streams <- parallel::clusterCall(cl[1], make_streams,
        initial_seed, length(stochastic))[[1]]
      for (j in seq_along(stochastic)) jobs[[stochastic[j]]]$seed <- seed_streams[[j]]
    }
    # Export once per worker; subsequent tasks carry only starts, methods and seeds.
    # This environment is local. clusterExport writes into the workers' workspaces.
    worker_data <- list2env(list(.rcaiman_cie_worker = list(fn = fcost, runner = runner)),
      parent = emptyenv())
    parallel::clusterExport(cl, ".rcaiman_cie_worker", envir = worker_data)
    dispatch <- function(job) {
      worker <- get(".rcaiman_cie_worker", envir = .GlobalEnv)
      worker$runner(job, worker$fn)
    }
    environment(dispatch) <- baseenv()
    fits <- parallel::parLapplyLB(cl, jobs, dispatch, chunk.size = 1L)
  } else {
    fits <- lapply(jobs, runner, fn = fcost)
  }
  # Assemble in declared order, irrespective of worker completion order.
  for (i in seq_along(start_vectors)) {
    start <- start_vectors[[i]]
    initial <- evaluate(start, start, initial_source)
    record(initial, i, initial_source, NA_character_, NA_integer_, "", fcost(start))
    for (j in seq_along(method)) {
      outcome <- fits[[(i-1L) * length(method) + j]]
      fit <- outcome$fit
      messages <- outcome$messages
      candidate <- if (!is.null(fit)) evaluate(fit$par, start, "optimized", method[j], fit) else NULL
      if (!is.null(fit) && is.null(candidate)) messages <- c(messages, "Unusable predictions or zenith reference.")
      record(candidate, i, "optimized", method[j],
        if (is.null(fit)) NA_integer_ else fit$convergence,
        paste(unique(messages), collapse = "; "), if (is.null(fit)) Inf else fit$value)
    }
  }
  attempts <- do.call(rbind, attempts)
  if (!any(is.finite(attempts$metric))) stop("No initial or optimized model produced a usable result.", call. = FALSE)
  best <- which.min(attempts$metric)
  result <- models[[best]]
  attempts$selected <- seq_len(nrow(attempts)) == best
  result$attempts <- attempts
  if (!any(is.finite(attempts$metric[attempts$source == "optimized"]))) {
    warning("No optimization attempt produced a usable result; returning the best initial model. See `attempts`.", call. = FALSE)
  } else if (result$source == "optimized" && result$opt_result$convergence != 0) {
    warning("The selected optimization did not converge. See `opt_result` and `attempts`.", call. = FALSE)
  }
  result
}

.cie_cached_predictor <- function(zenith, azimuth, sun_angles) {
  chi <- calc_spherical_distance(zenith, azimuth, sun_angles$zenith, sun_angles$azimuth)
  cos_zenith <- cos(zenith)
  cos_chi_squared <- cos(chi)^2
  sun_zenith <- sun_angles$zenith
  cos_sun_squared <- cos(sun_zenith)^2
  predict <- function(params) {
    gradation <- (1 + params[1] * exp(params[2] / cos_zenith)) /
      (1 + params[1] * exp(params[2]))
    indicatrix <- (1 + params[3] * (exp(params[4] * chi) - exp(params[4] * pi/2)) +
      params[5] * cos_chi_squared) /
      (1 + params[3] * (exp(params[4] * sun_zenith) - exp(params[4] * pi/2)) +
        params[5] * cos_sun_squared)
    unname(gradation * indicatrix)
  }
  environment(predict) <- list2env(list(chi = chi, cos_zenith = cos_zenith,
    cos_chi_squared = cos_chi_squared, sun_zenith = sun_zenith,
    cos_sun_squared = cos_sun_squared), parent = baseenv())
  predict
}

.cie_run_optim <- function(job, fn) {
  # Non-NULL seeds are used only inside dedicated worker processes.
  if (!is.null(job$seed)) assign(".Random.seed", job$seed, envir = .GlobalEnv)
  messages <- character()
  fit <- tryCatch(withCallingHandlers(
    stats::optim(par = job$start, fn = fn, method = job$method),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      messages <<- c(messages, conditionMessage(e))
      NULL
    })
  if (!is.null(fit$message)) messages <- c(messages, fit$message)
  list(fit = fit, messages = messages)
}

# Shared by fitting and the final refit in solar-position optimization.
.check_cie_general_sky_type <- function(general_sky_type) {
  .assert_choice(general_sky_type, c("Overcast", "Partly cloudy", "Clear"),
                 multiple = TRUE, allow_null = TRUE)
  if (!is.null(general_sky_type) && !length(general_sky_type) %in% 1:2) {
    stop("`general_sky_type` must contain one or two types; use NULL for all standard skies.",
         call. = FALSE)
  }
  invisible(TRUE)
}
