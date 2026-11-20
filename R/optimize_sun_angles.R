#' Optimize sun angular coordinates
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Refine a solar direction by repeatedly fitting the CIE sky coefficients at
#' candidate positions. Use the relative-radiance samples of an existing model.
#'
#' @details
#' A change in solar position generally requires a change in sky coefficients.
#' Each candidate is therefore evaluated with [fit_cie_model()], starting from
#' the coefficients of the current model. This is a search for a solar direction
#' supported by the sampled sky radiance, rather than a detection of the disk
#' in the image. The distribution and purity of sky samples affect the result.
#'
#' The solar search combines discrete exploration and continuous refinement
#' on the same spherical geometry. Exploration compares a finite set of
#' positions distributed around the current estimate. Refinement then lets
#' [stats::optim()] propose positions between and beyond those candidates.
#' Both stages refit the sky coefficients at every new position. A final
#' comparison broadens the choice of sky models at the resulting solar position.
#'
#' To describe a trial position, the search uses an *offset*: an angular
#' displacement from a reference solar direction, called the search center.
#' Imagine a plane touching the unit sphere at that direction. Two signed
#' components, in degrees along perpendicular directions in this tangent
#' plane, specify the displacement. Their combined length gives the angular
#' distance to travel on the sphere, and their orientation gives the direction
#' of travel along a great circle. Thus, offset `(0, 0)` is the center, and
#' offset `(3, 4)` reaches a position 5 deg away because
#' \eqn{\sqrt{3^2 + 4^2} = 5}.
#'
#' These components are not increments added to solar azimuth and zenith.
#' They are converted together into a new direction on the sphere. This
#' representation accommodates travel across azimuth zero and near the zenith,
#' where a change in azimuth alone is a poor measure of angular displacement.
#' Offsets are internal search variables; input and output solar positions
#' remain tables of `zenith` and `azimuth` in radians.
#'
#' 1. **Discrete exploration.** Generate Fibonacci candidates directly within a
#'    spherical cap around the current sun. Include its current position, and
#'    shrink the cap from `start_radius` to `tol_angle` over `max_iter`
#'    iterations. This stage chooses a finite set of offsets within each cap.
#'    Refit the sky with `method_sky` at each position. All candidates
#'    within an iteration start from the same coefficients. Accept an improvement
#'    and use it as the center of the next cap. A stationary estimate alone does
#'    not establish convergence: stop at the final cap radius, after
#'    `max_no_improve` consecutive iterations without improvement, or at
#'    `max_iter`. With `max_iter = 1`, evaluate the initial cap only.
#'
#' 2. **Continuous local refinement.** Starting at offset `(0, 0)`, each method
#'    in `method_sun` varies the two offset components through [stats::optim()].
#'    The components are continuous variables: trial positions need not belong
#'    to the Fibonacci pattern or remain within the last exploration cap.
#'    The same geometric conversion maps each offset to a solar direction,
#'    and another fit with `method_sky` supplies its cost. The optimizer
#'    therefore evaluates a position together with its refitted sky, rather
#'    than moving the sun under fixed sky coefficients. Offsets may reach
#'    at most 180 deg from the center; the solar domain and displacement
#'    penalty described below also apply. Repeated positions within a
#'    refinement reuse their evaluations. Each method has a budget of
#'    `max_sun_evaluations` distinct, admissible positions. Retain the best
#'    evaluated candidate, even if the optimizer fails or exhausts its budget.
#'
#' 3. **Final sky comparison.** At the resulting solar position, fit the 15
#'    standard skies again, or the starts selected by `custom_sky_coef`,
#'    `std_sky_no`, or `general_sky_type`. These selectors and `method` apply
#'    only to this stage. Retain the preceding model if the comparison does not
#'    improve its mean squared deviation. This wider comparison is performed
#'    once, rather than at every trial solar position.
#'
#' During exploration and local refinement, rank candidates by
#' \deqn{\frac{\sqrt{\mathrm{MSD}}}{\overline{L_r}}
#'        \max\left(1, \frac{\delta}{10\ \mathrm{deg}}\right),}
#' where MSD is the refitted model's mean squared deviation, \eqn{\overline{L_r}}
#' is mean observed relative radiance, and \eqn{\delta} is angular distance from
#' the current search center. This is a soft penalty: displacements beyond
#' 10 deg are allowed when their improvement compensates for the distance.
#' The center moves after each accepted exploration step or refinement method;
#' it is not tied to the initial solar detection. Within one refinement, both
#' the center and the initial coefficients remain fixed, so evaluation order
#' does not change the cost assigned to a position.
#'
#' The zenith DN reference is held fixed throughout the solar search.
#' `optim_zenith_dn` applies only to the final sky comparison; see
#' [fit_cie_model()] for its interpretation and effect on the comparison scale.
#' The coefficient fits use that function's existing physical penalties and
#' selection rule, including the horizon reference for the brightness penalty
#' when the sun is below the horizon. No additional constraints on sky
#' coefficients are introduced.
#'
#' Search positions are restricted to solar zenith angles from 0 to 96 deg,
#' including the civil-twilight range. This is the search domain, not a claim
#' that photographs beyond it are impossible. An initial position outside the
#' domain is rejected rather than silently changed. Candidates outside it are
#' excluded, so a cap near the boundary can contain fewer than `n_starts` rows.
#'
#' @param model list returned by [fit_cie_model()]. Solar and sample directions
#'   use the external reference frame and are in radians.
#' @param method_sun character vector or `NULL`. Optimization methods for continuous
#'   solar refinement using the two offset components explained in *Details*.
#'   Defaults to `"Nelder-Mead"`, to limit the number of nested sky fits. Supports
#'   the multivariate methods of [stats::optim()]. Use `NULL` to skip
#'   local refinement. Methods run in order, each from the best result so far.
#' @param max_iter nonnegative integer-like numeric vector of length one.
#'   Maximum number of exploration iterations; zero skips exploration.
#' @param tol_angle positive numeric vector of length one. Final exploration
#'   radius in degrees, no greater than `start_radius`.
#' @param max_no_improve positive integer-like numeric vector of length one.
#'   Maximum consecutive exploration iterations without improvement.
#' @param n_starts positive integer-like numeric vector of length one. Maximum
#'   positions per exploration iteration, including the current sun.
#' @param start_radius positive numeric vector of length one, at most 180.
#'   Initial exploration radius in degrees.
#' @param method_sky character vector. Coefficient-fitting methods used at every
#'   trial solar position, starting from the current coefficients. Defaults to
#'   `"Nelder-Mead"`. Supports `"Nelder-Mead"`, `"BFGS"`, `"CG"`, and
#'   `"L-BFGS-B"`. SANN is excluded here to keep repeated cost evaluations
#'   deterministic. Independent of the winning `model$method`.
#' @param max_sun_evaluations positive integer-like numeric vector of length one.
#'   Maximum distinct admissible solar positions fitted per local refinement
#'   method, including its center. Defaults to 60. Exploration and the final
#'   sky comparison are outside this budget.
#' @param method character vector. Coefficient-fitting methods for the final
#'   sky comparison, as in [fit_cie_model()].
#' @param optim_zenith_dn logical vector of length one. Optimize the zenith DN
#'   reference only in the final sky comparison; keep it fixed during the solar
#'   search. See [fit_cie_model()] for the reference-adjustment penalty.
#' @param parallel,cores,logical,leave_free Controls passed to [fit_cie_model()]
#'   for the final sky comparison only. The nested solar search is sequential.
#' @inheritParams fit_cie_model
#'
#' @return A model with the same components as [fit_cie_model()], plus:
#' \describe{
#'   \item{`method_sun`}{The last local-refinement method that improved the solar
#'     search, or `NULL` if none did.}
#'   \item{`sun_search`}{A list with `history`, `stop_reason`, and `refinements`.
#'     `history` records evaluated positions (`zenith` and `azimuth` in radians),
#'     stage, iteration, method, MSD (`metric`), penalized cost (`score`),
#'     acceptance, and diagnostic messages. `stop_reason` describes the end of
#'     exploration. `refinements` records each local optimizer's convergence,
#'     number of evaluated positions, budget exhaustion, and messages.
#'     Failed candidate fits have infinite cost. `accepted` marks successive
#'     retained models, not only the final result.}
#' }
#'
#' @note
#' The search is local and does not guarantee a global optimum. Some skies
#' provide little information about solar position; for example, a CIE model
#' with both indicatrix coefficients `c` and `e` equal to zero is independent
#' of the sun. Equal-cost candidates retain the existing position.
#'
#' Candidate generation and `method_sky` are deterministic. SANN in
#' `method_sun` or the final `method` uses randomness; set a seed in the calling
#' script for repeatability. This function does not set a seed or change the
#' caller's random-number generator configuration. Finite nonconverged sky
#' fits may be used, as in [fit_cie_model()]; inspect the search diagnostics.
#'
#' @seealso [fit_cie_model()], [calc_relative_radiance()],
#'   [optimize_sky_mask()], [optimize_sampling_distance()]
#' @export
#'
#' @examples
#' # Known synthetic sky, with a displaced initial solar position.
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(20, 15),
#'   normalization_radius = 18,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(40, 30), geometry)
#' d <- direction_image(g)
#' sun <- data.frame(zenith = 45 * pi/180, azimuth = 30 * pi/180)
#' sky <- cie_image(d, sun, as.numeric(cie_table[12, 1:5])) * 1000
#' names(sky) <- "Blue"
#' indices <- expand.grid(col = seq(2, 40, 4), row = seq(2, 30, 4))
#' indices <- rem_unselected(indices, !is.na(sky))
#' samples <- add_directions(extract_dn(sky, indices, use_window = FALSE), d = d)
#' rr <- calc_relative_radiance(samples, "Blue", zenith_dn = 1000)
#' initial_sun <- sun
#' initial_sun$zenith <- 48 * pi/180
#' fit <- fit_cie_model(rr, initial_sun, std_sky_no = 12, method = "Nelder-Mead")
#' refined <- optimize_sun_angles(fit, max_iter = 2, n_starts = 5,
#'   max_sun_evaluations = 12, std_sky_no = 12, method = "Nelder-Mead")
#' rbind(initial = fit$sun_angles, refined = refined$sun_angles, known = sun)
#' c(initial = fit$metric, refined = refined$metric)
#' refined$sun_search$history
#'
#' \dontrun{
#' # Broader final comparison, with all 15 skies and the default fit methods.
#' set.seed(7)
#' refined <- optimize_sun_angles(fit, parallel = TRUE, cores = 4)
#' fitted_sky <- cie_image(d, refined$sun_angles, refined$coef) *
#'   refined$relative_radiance$zenith_dn
#' plot(fitted_sky)
#' }
optimize_sun_angles <- function(model,
                                method_sun = "Nelder-Mead",
                                max_iter = 20,
                                tol_angle = 0.1,
                                max_no_improve = 5,
                                n_starts = 10,
                                start_radius = 5,
                                custom_sky_coef = NULL,
                                std_sky_no = NULL,
                                general_sky_type = NULL,
                                optim_zenith_dn = FALSE,
                                method_sky = "Nelder-Mead",
                                max_sun_evaluations = 60,
                                method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
                                parallel = FALSE,
                                cores = NULL,
                                logical = TRUE,
                                leave_free = 1) {

  .assert_cie_model(model)
  .check_cie_sun_angles(model$sun_angles)
  if (model$sun_angles$zenith > 96 * pi/180) {
    stop("The initial solar zenith angle must not exceed 96 degrees.", call. = FALSE)
  }
  .check_vector(model$coef, "numeric", 5)
  .check_vector(model$metric, "numeric", 1, sign = "nonnegative")
  .check_vector(max_iter, "integerish", 1, sign = "nonnegative")
  for (value in c("max_no_improve", "n_starts", "max_sun_evaluations")) {
    .check_vector(get(value), "integerish", 1, sign = "positive", name = value)
  }
  .check_vector(tol_angle, "numeric", 1, sign = "positive")
  .check_vector(start_radius, "numeric", 1, sign = "positive")
  if (start_radius > 180 || tol_angle > start_radius) {
    stop("Require `tol_angle <= start_radius <= 180` degrees.", call. = FALSE)
  }
  supported <- c("Nelder-Mead", "BFGS", "CG", "SANN", "L-BFGS-B")
  for (value in c("method_sun", "method_sky", "method")) {
    x <- get(value)
    if (value == "method_sun" && is.null(x)) next
    allowed <- if (value == "method_sky") setdiff(supported, "SANN") else supported
    if (!is.character(x) || anyNA(x) || any(!x %in% allowed) ||
        (value != "method_sun" && !length(x))) {
      stop("`", value, "` must contain supported optimization methods.", call. = FALSE)
    }
  }
  method_sun <- unique(method_sun)
  method_sky <- unique(method_sky)
  .check_vector(optim_zenith_dn, "logical", 1)
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, sign = "positive", allow_null = TRUE)
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")
  .check_cie_general_sky_type(general_sky_type)
  if (!is.null(custom_sky_coef)) {
    if (!is.null(std_sky_no) || !is.null(general_sky_type)) {
      stop("Do not combine `custom_sky_coef` with standard sky selectors.", call. = FALSE)
    }
    x <- custom_sky_coef
    if (is.numeric(x) && is.null(dim(x)) && length(x) == 5L) x <- matrix(x, nrow = 1)
    if (!is.matrix(x) || !is.numeric(x) || ncol(x) != 5L || !nrow(x) || any(!is.finite(x))) {
      stop("`custom_sky_coef` must contain finite five-coefficient starts.", call. = FALSE)
    }
  }
  if (!is.null(std_sky_no)) {
    .check_vector(std_sky_no, "integerish", sign = "positive")
    if (!length(std_sky_no) || any(std_sky_no > 15)) {
      stop("`std_sky_no` must contain numbers from 1 to 15.", call. = FALSE)
    }
  }
  rr <- model$relative_radiance
  .check_vector(rr$zenith_dn, "numeric", 1, sign = "positive")
  .check_columns(rr$sky_points, c("zenith", "azimuth", "relative_radiance"))
  for (column in c("zenith", "azimuth", "relative_radiance")) {
    .check_vector(rr$sky_points[[column]], "numeric", name = column)
  }
  if (!nrow(rr$sky_points) || any(rr$sky_points$zenith < 0 | rr$sky_points$zenith > pi/2) ||
      any(rr$sky_points$azimuth < 0 | rr$sky_points$azimuth >= 2*pi) ||
      !is.finite(mean(rr$sky_points$relative_radiance)) || mean(rr$sky_points$relative_radiance) <= 0) {
    stop("Supply finite sky samples with external directions in radians and positive mean relative radiance.", call. = FALSE)
  }
  mean_radiance <- mean(rr$sky_points$relative_radiance)
  history <- list()
  record <- function(sun, fit, stage, iteration, method, score, message = "") {
    id <- length(history) + 1L
    history[[id]] <<- data.frame(stage = stage, iteration = iteration,
      method = method, zenith = sun$zenith, azimuth = sun$azimuth,
      metric = if (is.null(fit)) Inf else fit$metric, score = score,
      accepted = FALSE, message = message, stringsAsFactors = FALSE)
    id
  }
  accept <- function(id) history[[id]]$accepted <<- TRUE
  score <- function(fit, center) {
    delta <- calc_spherical_distance(fit$sun_angles$zenith, fit$sun_angles$azimuth,
      center$zenith, center$azimuth) * 180/pi
    sqrt(fit$metric) / mean_radiance * max(1, delta/10)
  }
  attempt <- function(sun, anchor, final = FALSE) {
    messages <- character()
    fit <- tryCatch(withCallingHandlers({
      if (final) {
        fit_cie_model(anchor$relative_radiance, sun,
          custom_sky_coef = custom_sky_coef, std_sky_no = std_sky_no,
          general_sky_type = general_sky_type, optim_zenith_dn = optim_zenith_dn,
          method = method, parallel = parallel, cores = cores,
          logical = logical, leave_free = leave_free)
      } else {
        fit_cie_model(anchor$relative_radiance, sun,
          custom_sky_coef = anchor$coef, method = method_sky)
      }
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
      messages <<- c(messages, conditionMessage(e)); NULL
    })
    list(fit = fit, message = paste(unique(messages), collapse = "; "))
  }
  id <- record(model$sun_angles, model, "initial", 0L, NA_character_,
    sqrt(model$metric)/mean_radiance)
  accept(id)
  stop_reason <- if (max_iter == 0) "exploration_skipped" else "max_iter"
  no_improve <- 0L
  selected_method <- NULL
  radii <- if (max_iter > 1) exp(seq(log(start_radius), log(tol_angle), length.out = max_iter)) else start_radius
  for (iteration in seq_len(max_iter)) {
    anchor <- model
    candidates <- .sun_cap_candidates(anchor$sun_angles, n_starts, radii[iteration])
    best <- model; best_score <- sqrt(model$metric)/mean_radiance; best_id <- NULL
    for (j in seq_len(nrow(candidates))) {
      sun <- candidates[j, , drop = FALSE]
      trial <- attempt(sun, anchor)
      cost <- if (is.null(trial$fit)) Inf else score(trial$fit, anchor$sun_angles)
      id <- record(sun, trial$fit, "exploration", iteration, NA_character_, cost, trial$message)
      if (cost < best_score) {
        best <- trial$fit; best_score <- cost; best_id <- id
      }
    }
    if (!is.null(best_id)) {
      model <- best; accept(best_id); no_improve <- 0L
    } else no_improve <- no_improve + 1L
    if (radii[iteration] <= tol_angle * (1 + 1e-12)) {
      stop_reason <- "tol_angle"; break
    }
    if (no_improve >= max_no_improve) {
      stop_reason <- "max_no_improve"; break
    }
  }
  refinements <- list()
  for (solar_method in method_sun) {
    anchor <- model
    best <- model; best_score <- sqrt(model$metric)/mean_radiance; best_id <- NULL
    cache <- new.env(parent = emptyenv())
    evaluations <- 0L; exhausted <- FALSE
    objective <- function(offset) {
      if (any(!is.finite(offset)) || sqrt(sum(offset^2)) > 180) return(Inf)
      sun <- .sun_offset_direction(anchor$sun_angles, offset[1], offset[2])
      if (sun$zenith > 96*pi/180) return(Inf)
      key <- paste(sprintf("%.17g", unlist(sun)), collapse = ":")
      if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
      if (evaluations >= max_sun_evaluations) {
        exhausted <<- TRUE
        return(Inf)
      }
      evaluations <<- evaluations + 1L
      trial <- attempt(sun, anchor)
      cost <- if (is.null(trial$fit)) Inf else score(trial$fit, anchor$sun_angles)
      id <- record(sun, trial$fit, "refinement", 0L, solar_method, cost, trial$message)
      assign(key, cost, cache)
      if (cost < best_score) {
        best <<- trial$fit; best_score <<- cost; best_id <<- id
      }
      cost
    }
    messages <- character()
    opt <- tryCatch(withCallingHandlers(
      stats::optim(c(0, 0), objective, method = solar_method,
        control = list(maxit = max_sun_evaluations, parscale = rep(start_radius, 2))),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) {
        messages <<- c(messages, conditionMessage(e)); NULL
      })
    if (!is.null(opt$message)) messages <- c(messages, opt$message)
    if (!is.null(best_id)) {
      model <- best; accept(best_id); selected_method <- solar_method
    }
    refinements[[length(refinements) + 1L]] <- data.frame(method = solar_method,
      convergence = if (is.null(opt)) NA_integer_ else opt$convergence,
      evaluations = evaluations, budget_exhausted = exhausted,
      message = paste(unique(messages), collapse = "; "), stringsAsFactors = FALSE)
  }
  trial <- attempt(model$sun_angles, model, final = TRUE)
  final_score <- if (is.null(trial$fit)) Inf else sqrt(trial$fit$metric)/mean_radiance
  id <- record(model$sun_angles, trial$fit, "final", 0L, NA_character_, final_score, trial$message)
  if (!is.null(trial$fit) && trial$fit$metric < model$metric) {
    model <- trial$fit; accept(id)
  }
  if (is.null(trial$fit)) {
    warning("The final sky comparison failed; retaining the preceding model. See `sun_search$history`.", call. = FALSE)
  }
  model["method_sun"] <- list(selected_method)
  model$sun_search <- list(history = do.call(rbind, history), stop_reason = stop_reason,
    refinements = if (length(refinements)) do.call(rbind, refinements) else
      data.frame(method = character(), convergence = integer(), evaluations = integer(),
        budget_exhausted = logical(), message = character()))
  model
}

# Exponential map from a two-dimensional angular offset (degrees).
.sun_offset_direction <- function(center, east, south) {
  theta <- center$zenith; phi <- center$azimuth
  v <- c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
  e <- c(-sin(phi), cos(phi), 0)
  s <- c(cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta))
  rho <- sqrt(east^2 + south^2) * pi/180
  if (rho == 0) return(center[, c("zenith", "azimuth"), drop = FALSE])
  xyz <- cos(rho)*v + sin(rho)/rho * (east*e + south*s)*pi/180
  azimuth <- atan2(xyz[2], xyz[1]) %% (2*pi)
  if (azimuth >= 2*pi) azimuth <- 0
  data.frame(zenith = atan2(sqrt(sum(xyz[1:2]^2)), xyz[3]),
    azimuth = azimuth)
}

.sun_cap_candidates <- function(center, n, radius) {
  candidates <- center[, c("zenith", "azimuth"), drop = FALSE]
  if (n == 1) return(candidates)
  j <- seq_len(n - 1L)
  # Area-uniform radii, computed without cancellation in very small caps.
  rho <- 2*asin(sqrt((j-.5)/(n-1L))*sin(radius*pi/360)) * 180/pi
  phi <- 2*pi*(j-1)/((1+sqrt(5))/2)
  extra <- lapply(j, function(i) .sun_offset_direction(center,
    rho[i]*cos(phi[i]), rho[i]*sin(phi[i])))
  candidates <- rbind(candidates, do.call(rbind, extra))
  candidates <- candidates[candidates$zenith <= 96*pi/180, , drop = FALSE]
  rownames(candidates) <- NULL
  candidates
}
