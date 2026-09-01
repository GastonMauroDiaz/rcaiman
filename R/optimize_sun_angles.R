# .plot_it <- function(starts) {
#   .normalize_angles <- function(zen_deg, azi_deg) {
#     if (zen_deg < 0) {
#       zen_deg <- -zen_deg      # Make zenith angle positive
#       azi_deg <- azi_deg + 180 # Rotate azimuth angle 180 deg
#     }
#
#     # Ensure azimuth angle is in [0, 360)
#     azi_deg <- azi_deg %% 360
#
#     angles <- c(zen_deg, azi_deg)
#     names(angles) <- c("zen_deg", "azi_deg")
#     angles
#   }
#   za <- lapply(1:nrow(starts), function(i)
#     .normalize_angles(starts[i,1], starts[i,2]))
#   zen <- unlist(lapply(za, function(x) x[1])) * pi/180
#   azi <- unlist(lapply(za, function(x) x[2])) * pi/180
#   row_col <- row_col_from_zen_azi(z, a, zen, azi)
#   plot(row_col, asp = 1)
# }

.compute_spacing <- function(n, radius_deg) {
  r <- radius_deg * pi / 180

  frac <- 1 - cos(r)
  frac <- max(frac, 1e-6)

  N <- n / frac

  d <- sqrt(2 * pi / N)
  theta <- 2 * asin(pmin(1, d / 2))

  list(
    spacing = theta * 180 / pi,
    N_full = 4 * pi / (d^2)
  )
}

.generate_starts <- function(center, n, radius_deg, grid_xyz) {

  if (n == 1) return(matrix(center, nrow = 1))

  center_xyz <- .directions_to_cartesian(
    data.frame(zen = center["zen_deg"] %>% .degree2radian(),
               azi = center["azi_deg"] %>% .degree2radian())
  )

  nn <- .nn_spherical_fast(
    xyz = grid_xyz,
    xyz0 = c(x = center_xyz$x, y = center_xyz$y, z = center_xyz$z),
    k = n - 1,
    search_radius = .degree2radian(radius_deg)
  )

  pts <- grid_xyz[nn$index, ]

  ang <- .cartesian_to_directions(as.matrix(pts))

  starts <- rbind(center,
                  cbind(zen_deg = ang[, "zen"] %>% .radian2degree(),
                        azi_deg = ang[, "azi"] %>% .radian2degree()))

  starts
}

.optimize_sun_angles <- function(sun_angles, rr, sky_coef, method_sun) {
  .normalize_angles <- function(zen_deg, azi_deg) {
    if (zen_deg < 0) {
      zen_deg <- -zen_deg      # Make zenith angle positive
      azi_deg <- azi_deg + 180 # Rotate azimuth angle 180 deg
    }

    # Ensure azimuth angle is in [0, 360)
    azi_deg <- azi_deg %% 360

    angles <- c(zen_deg, azi_deg)
    names(angles) <- c("zen_deg", "azi_deg")
    angles
  }
  .refine_sun_coord <- function(param) {
    zen_deg <- param[1]
    azi_deg <- param[2]

    chi <- calc_spherical_distance(zen_deg %>% .degree2radian(),
                                   azi_deg %>% .degree2radian(),
                                   sun_angles["zen_deg"] %>% .degree2radian(),
                                   sun_angles["azi_deg"] %>% .degree2radian())


    pred <- .cie_sky_model(AzP = rr$sky_points$azi %>% .degree2radian(),
                           Zp = rr$sky_points$zen %>% .degree2radian(),
                           AzS = azi_deg %>% .degree2radian(),
                           Zs =  zen_deg %>% .degree2radian(),
                           sky_coef[1], sky_coef[2], sky_coef[3],
                           sky_coef[4], sky_coef[5])

    (.calc_rmse(pred - rr$sky_points$rr) / mean(rr$sky_points$rr)) *
      max(1, .radian2degree(chi)/10)
  }
  fit <- tryCatch(stats::optim(c(sun_angles["zen_deg"],
                                 sun_angles["azi_deg"]),
                               .refine_sun_coord,
                               method = method_sun),
                  error = function(e) NULL)

  tryCatch(.normalize_angles(fit$par[1], fit$par[2]),
           error = function(e) NULL)
}

#' Optimize sun angular coordinates
#'
#' @description
#' Refine sun zenith and azimuth angles to minimize mismatch between observed
#' and predicted relative radiance under a CIE sky model.
#'
#' @details
#' The procedure combines an iterative heuristic search with a final local
#' optimization step. The heuristic stage explores the parameter space to avoid
#' poor local minima, while the final step refines the best solution found.
#'
#' 1. **Heuristic candidate generation**
#'    A spherical Fibonacci lattice is generated at each iteration with angular
#'    spacing adapted to the current search radius. The radius decreases
#'    monotonically across iterations, while the sampling density increases
#'    accordingly. This defines a coarse-to-fine exploration of the parameter
#'    space with approximately constant sampling effort.
#'
#'    The angular domain extends beyond the horizon (up to ~96 deg zenith
#'    angle), allowing the method to operate under crepuscular illumination
#'    conditions. The upper bound (~96 deg) corresponds to civil twilight.
#'    Beyond this range, image acquisition is not feasible due to insufficient
#'    illumination, and therefore such conditions are not considered.
#'
#' 2. **Heuristic model evaluation**
#'    Each candidate position is evaluated by fitting a CIE sky model using
#'    [fit_cie_model()] with fixed coefficients as starting values
#'    (`custom_sky_coef = model$coef`). The candidate yielding the lowest model
#'    metric is selected.
#'
#' 3. **Iterative update**
#'    If an improved solution is found, the model is updated and used as the
#'    reference for the next iteration. This defines an iterative search process
#'    guided by local exploration of the angular neighborhood.
#'
#' 4. **Convergence criteria**
#'    Iteration stops when any of the following conditions is met:
#'    \itemize{
#'      \item angular change between successive estimates is below `tol_angle`,
#'      \item no improvement is observed for `max_no_improve` iterations,
#'      \item the number of iterations reaches `max_iter`.
#'    }
#'
#'    Iteration may also stop when the angular resolution becomes finer than
#'    required by `tol_angle`, or when the implied sampling density would lead
#'    to excessive computational cost.
#'
#' 5. **Final local refinement**
#'    After convergence, a local optimization is performed from the final
#'    estimate using each method in `method_sun` via [stats::optim()]. The
#'    resulting candidates are evaluated with [fit_cie_model()], and the best
#'    solution is retained.
#'
#' This strategy separates global exploration (heuristic stage) from local
#' exploitation (final refinement), improving robustness against local minima.
#'
#' @param model list returned by [fit_cie_model()].
#' @param method_sun character vector. Optimization methods supported by
#'   [stats::optim()]. Used only in the final refinement step.
#' @param max_iter integer of length one. Maximum number of outer iterations.
#' @param tol_angle numeric of length one. Angular tolerance in degrees.
#' @param max_no_improve integer of length one. Maximum number of consecutive
#'   iterations without improvement.
#' @param n_starts integer of length one. Number of candidate sun positions
#'   evaluated per iteration. Must be >= 1.
#' @param start_radius numeric of length one. Angular radius (deg) defining the
#'   neighborhood used to generate candidate positions.
#'
#' @inheritParams fit_cie_model
#'
#' @return List like `model` with updated components and an additional element:
#' \describe{
#'   \item{`method_sun`}{character string indicating the selected optimization
#'   method in the final refinement step, or `NULL` if no improvement was found.}
#' }
#'
#' @note
#' The heuristic stage is deterministic given `n_starts` and `start_radius`,
#' as it relies on a fixed spherical Fibonacci lattice and nearest-neighbor
#' selection in angular space.
#'
#' Sun angles are internally normalized to enforce `zen_deg >= 0` and
#' `azi_deg` in [0, 360).
#'
#' In contrast to [estimate_sun_angles()], this function allows zenith angles
#' greater than 90 deg to support crepuscular conditions.
#'
#' @export
#'
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
#' plot(r)
#'
#' # See fit_cie_model() for details on below file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' sun_angles <- c(zen_deg = 49.5, azi_deg = 33.4)
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#' rr <- calc_relative_radiance(sky_points, col_id = names(r))
#'
#' set.seed(7)
#' model <- fit_cie_model(rr, sun_angles, general_sky_type = "Clear")
#'
#' print(model$sun_angles)
#' print(model$metric)
#'
#' plot(model$rr$sky_points$rr, model$rr$sky_points$pred)
#' abline(0,1)
#' lm(model$rr$sky_points$pred~model$rr$sky_points$rr) %>% summary()
#'
#' model2 <- optimize_sun_angles(model,
#'                               max_iter = 30,
#'                               tol_angle = 0.1,
#'                               max_no_improve = 20,
#'                               n_starts = 30,
#'                               start_radius = 10
#' )
#'
#' plot(model2$rr$sky_points$rr, model2$rr$sky_points$pred)
#' abline(0,1)
#'
#' print(model$sun_angles)
#' print(model2$sun_angles)
#' print(model$metric)
#' print(model2$metric)
#' model2$method_sun
#'
#' plot(r)
#' .add_sun_disk <- function(sun_angles,  col) {
#'   sun_row_col <- row_col_from_zen_azi(z, a,
#'                                       sun_angles["zen_deg"] * pi/180,
#'                                       sun_angles["azi_deg"] * pi/180)
#'   points(sun_row_col[2], nrow(r) - sun_row_col[1],
#'          col = col, pch = 8, cex = 3)
#' }
#' .add_sun_disk(sun_angles, col = "yellow")
#' .add_sun_disk(model2$sun_angles, col = "orange")
#' }
optimize_sun_angles <- function(model,
                                method_sun = c("Nelder-Mead", "BFGS", "CG", "SANN"),
                                max_iter = 20,
                                tol_angle = 0.1,
                                max_no_improve = 5,
                                n_starts = 10,
                                start_radius = 5,
                                custom_sky_coef = NULL,
                                std_sky_no = NULL,
                                general_sky_type = NULL,
                                optim_zenith_dn = FALSE) {

  .assert_cie_model(model)

  if (model$sun_angles["zen_deg"] > 96) {
    message(
      sprintf(
        "Provided sun zenith angle was %f. Since it is unrealistic, it was changed to 96 deg (civic twilight) before optimization.",
        model$sun_angles["zen_deg"]
      )
    )
    model$sun_angles["zen_deg"] <- 96
  }

  my_col_names <- names(model$rr$sky_points)
  if ("dn" %in% my_col_names) {
    model$rr$sky_points <- model$rr$sky_points[, -grep("dn", my_col_names)]
  }

  iter <- 0
  no_improve <- 0
  best_metric <- model$metric

  # multi-scale schedule
  alpha <- (tol_angle / start_radius)^(1 / (max_iter - 1))

  # internal guardrail (not exposed)
  N_max_internal <- 2e5

  repeat {
    iter <- iter + 1
    if (iter > max_iter) break

    # radius decay
    radius_iter <- start_radius * alpha^(iter - 1)

    # adaptive spacing
    sp <- .compute_spacing(n_starts, radius_iter)
    spacing <- sp$spacing
    N_full <- sp$N_full

    # stopping criteria (resolution + cost)
    if (spacing < tol_angle / 2) break
    if (N_full > N_max_internal) break

    # regenerate grid at current scale
    grid_xyz <- .fibonacci(spacing, max_zenith_angle = 96) #civil twilight
    grid_xyz <- data.frame(grid_xyz)
    names(grid_xyz) <- c("x", "y", "z")

    pre_optim_sun <- model$sun_angles

    starts <- .generate_starts(
      model$sun_angles,
      n_starts,
      radius_iter,
      grid_xyz
    )

    models <- lapply(seq_len(nrow(starts)), function(i) {
      sun_angles <- starts[i, ]
      names(sun_angles) <- c("zen_deg", "azi_deg")

      fit_cie_model(
        model$rr,
        sun_angles,
        custom_sky_coef = model$coef,
        method = model$method
      )
    })

    metrics <- vapply(models, function(m) m$metric, numeric(1))
    i <- which.min(metrics)

    if (metrics[i] < model$metric) {
      model <- models[[i]]
    }

    delta_chi <- calc_spherical_distance(
      pre_optim_sun["zen_deg"],
      pre_optim_sun["azi_deg"],
      model$sun_angles["zen_deg"],
      model$sun_angles["azi_deg"]
    )

    if (model$metric < best_metric) {
      best_metric <- model$metric
      no_improve <- 0
    } else {
      no_improve <- no_improve + 1
    }

    if (delta_chi < tol_angle * pi / 180) break
    if (no_improve >= max_no_improve) break
  }

  suns <- lapply(method_sun, function(x) {
    .optimize_sun_angles(
      sun_angles = model$sun_angles,
      rr = model$rr,
      sky_coef = model$coef,
      method_sun = x
    )
  })

  valid <- vapply(suns, function(x) !is.null(x), logical(1))

  if (any(valid)) {
    method_sun_subset <- method_sun[valid]
    suns <- suns[valid]

    models <- lapply(suns, function(sun_angles)
      fit_cie_model(model$rr,
                    sun_angles,
                    custom_sky_coef = custom_sky_coef,
                    std_sky_no = std_sky_no,
                    general_sky_type = general_sky_type,
                    optim_zenith_dn = optim_zenith_dn,
                    method = model$method)
    )

    metrics <- vapply(models, function(m) m$metric, numeric(1))
    i <- which.min(metrics)

    if (metrics[i] < model$metric) {
      model <- models[[i]]
      .method_sun <- method_sun_subset[i]
    }
  }

  model$method_sun <- .method_sun
  model
}
