cie_validation_fixture <- function() {
  samples <- expand.grid(zenith = seq(.1, 1.4, length.out = 5),
                         azimuth = seq(0, 5.5, length.out = 6))
  sun <- data.frame(zenith = pi/4, azimuth = .3)
  coef <- as.numeric(cie_table[4, 1:5])
  samples$relative_radiance <- .cie_sky_model(samples$azimuth, samples$zenith,
    sun$azimuth, sun$zenith, coef[1], coef[2], coef[3], coef[4], coef[5])
  samples$id <- seq_len(nrow(samples))
  samples$pred <- -999
  list(relative_radiance = list(zenith_dn = 1200, sky_points = samples),
    sun_angles = sun, coef = coef, opt_result = NULL, method = NA_character_,
    start = coef, metric = 0)
}

validation_predictions <- function(model, k = 3, fold_id = NULL, ...) {
  args <- list(model = model, k = k, fold_id = fold_id,
    custom_sky_coef = NULL, std_sky_no = 4, general_sky_type = NULL,
    method = "Nelder-Mead", parallel = FALSE, cores = NULL,
    logical = TRUE, leave_free = 1)
  args[names(list(...))] <- list(...)
  do.call(.cie_validation_predictions, args)
}

test_that("held-out predictions use each training model and retain row order", {
  model <- cie_validation_fixture()
  original <- model
  calls <- list()
  coefs <- lapply(c(3, 5, 7), function(i) as.numeric(cie_table[i, 1:5]))
  local_mocked_bindings(fit_cie_model = function(relative_radiance, sun_angles,
    custom_sky_coef, std_sky_no, general_sky_type, optim_zenith_dn, method,
    parallel, cores, logical, leave_free) {
    j <- length(calls) + 1L
    calls[[j]] <<- relative_radiance
    expect_equal(sun_angles, model$sun_angles)
    expect_equal(relative_radiance$zenith_dn, 1200)
    expect_false(optim_zenith_dn)
    expect_null(custom_sky_coef)
    expect_equal(std_sky_no, 4)
    expect_equal(method, "Nelder-Mead")
    expect_false("pred" %in% names(relative_radiance$sky_points))
    list(coef = coefs[[j]])
  })
  folds <- rep(c(20, 8, 14), length.out = 30)
  result <- validation_predictions(model, k = NA, fold_id = folds)
  expected <- numeric(30)
  p <- model$relative_radiance$sky_points
  for (j in seq_along(calls)) {
    ids <- which(folds == unique(folds)[j])
    expect_equal(calls[[j]]$sky_points$id, setdiff(p$id, ids))
    a <- coefs[[j]]
    expected[ids] <- .cie_sky_model(p$azimuth[ids], p$zenith[ids],
      model$sun_angles$azimuth, model$sun_angles$zenith,
      a[1], a[2], a[3], a[4], a[5])
  }
  expect_equal(result$pred, expected)
  expect_equal(result$obs, p$relative_radiance)
  expect_equal(result$fold_id, folds)
  expect_named(result$fold_models, c("20", "8", "14"))
  expect_identical(model, original)
})

test_that("real fold fits recover a standard sky without indices or random starts", {
  model <- cie_validation_fixture()
  set.seed(910)
  seed <- .Random.seed
  result <- validation_predictions(model)
  expect_identical(.Random.seed, seed)
  expect_equal(result$pred, result$obs, tolerance = 1e-12)
  expect_identical(result$fold_id, rep(1:3, length.out = 30))
  expect_true(all(vapply(result$fold_models,
    function(x) x$relative_radiance$zenith_dn == 1200, logical(1))))
  expect_true(all(vapply(result$fold_models,
    function(x) identical(x$source, "standard"), logical(1))))
  # The coefficients and selected method of the original model are not starts.
  model$coef <- rep(99, 5)
  model$method <- "not an optimization method"
  expect_equal(validation_predictions(model)$pred, result$pred)
})

test_that("invalid folds and directions fail explicitly", {
  model <- cie_validation_fixture()
  expect_error(validation_predictions(model, k = 1), "`k`")
  expect_error(validation_predictions(model, k = 31), "`k`")
  expect_error(validation_predictions(model, k = 2.5), "`k`")
  expect_error(validation_predictions(model, fold_id = rep(1, 30)), "two folds")
  expect_error(validation_predictions(model, fold_id = 1:2), "fold_id")
  expect_error(validation_predictions(model, fold_id = c(NA, rep(1, 29))), "fold_id")
  expect_error(validation_predictions(model, fold_id = c(0, rep(1, 29))), "fold_id")
  model$relative_radiance$sky_points$zenith[1] <- pi
  expect_error(validation_predictions(model), "radians")
  model$relative_radiance$sky_points <- model$relative_radiance$sky_points[1, ]
  expect_error(validation_predictions(model), "two sky samples")
})

test_that("failures identify the fold rather than reuse the complete model", {
  model <- cie_validation_fixture()
  local_mocked_bindings(fit_cie_model = function(...) stop("no usable model"))
  expect_error(validation_predictions(model), "Fold 1: no usable model")
})

test_that("non-finite held-out predictions are not silently omitted", {
  model <- cie_validation_fixture()
  local_mocked_bindings(fit_cie_model = function(...) list(coef = rep(Inf, 5)))
  expect_error(validation_predictions(model), "Fold 1 produced non-finite")
})

test_that("public summaries use observed-on-predicted regression and retain all rows", {
  pred <- seq_len(12)
  obs <- pred + c(-.5, .2, -.1, .3, -.2, .1, -.3, .4, -.4, .5, 0, 20)
  local_mocked_bindings(.cie_validation_predictions = function(...) {
    list(pred = pred, obs = obs, fold_id = rep(1:3, 4), fold_models = list())
  })
  result <- validate_cie_model(NULL)
  expect_identical(result$pred, pred)
  expect_identical(result$obs, obs)
  expect_identical(result$is_outlier, c(rep(FALSE, 11), TRUE))
  keep <- !result$is_outlier
  expect_equal(unname(coef(result$lm)), unname(coef(lm(obs[keep] ~ pred[keep]))))
  expect_equal(result$metric, mean((obs[keep] - pred[keep])^2))
  expect_equal(sum(result$components), result$metric)
  expect_equal(result$rmse^2, result$metric)
  expect_equal(result$mae, mean(abs(obs[keep] - pred[keep])))
  expect_equal(result$medae, median(abs(obs[keep] - pred[keep])))
  expect_equal(result$r_squared, cor(pred[keep], obs[keep])^2)
})

test_that("laxity controls evaluation without changing held-out predictions", {
  pred <- seq_len(6)
  obs <- pred + c(0, 1, 2, 3, 4, 8)
  local_mocked_bindings(.cie_validation_predictions = function(...) {
    list(pred = pred, obs = obs, fold_id = rep(1:3, 2), fold_models = list())
  })
  default <- validate_cie_model(NULL)
  stricter <- validate_cie_model(NULL, laxity = 2)
  unfiltered <- validate_cie_model(NULL, laxity = NULL)
  expect_false(any(default$is_outlier))
  expect_identical(stricter$is_outlier, c(rep(FALSE, 5), TRUE))
  expect_identical(stricter$pred, default$pred)
  expect_identical(stricter$obs, default$obs)
  expect_identical(stricter$fold_id, default$fold_id)
  expect_identical(stricter$fold_models, default$fold_models)
  expect_equal(stricter$metric, mean((obs - pred)[1:5]^2))
  expect_false(any(unfiltered$is_outlier))
  expect_equal(unfiltered$metric, mean((obs - pred)^2))
})

test_that("invalid laxity fails before fitting folds", {
  local_mocked_bindings(.cie_validation_predictions = function(...) stop("Should not fit"))
  for (bad in list(0, -1, NA_real_, Inf, c(2, 3), "3")) {
    expect_error(validate_cie_model(NULL, laxity = bad), "laxity")
  }
})

test_that("zero MAD and constant observations remain evaluable", {
  local_mocked_bindings(.cie_validation_predictions = function(...) {
    list(pred = rep(1, 6), obs = rep(2, 6),
      fold_id = rep(1:3, 2), fold_models = list())
  })
  result <- validate_cie_model(NULL)
  expect_false(any(result$is_outlier))
  expect_equal(result$rmse, 1)
  expect_equal(result$mae, 1)
  expect_equal(result$medae, 1)
  expect_true(is.na(result$r_squared))
  expect_equal(result$components, c(SB = 1, NU = 0, LC = 0))
})

test_that("explicit fitting and parallel controls reach every training fit", {
  model <- cie_validation_fixture()
  calls <- 0L
  local_mocked_bindings(fit_cie_model = function(relative_radiance, sun_angles,
    custom_sky_coef, std_sky_no, general_sky_type, optim_zenith_dn, method,
    parallel, cores, logical, leave_free) {
    calls <<- calls + 1L
    expect_equal(custom_sky_coef, model$coef)
    expect_null(std_sky_no)
    expect_null(general_sky_type)
    expect_equal(method, c("BFGS", "CG"))
    expect_true(parallel)
    expect_equal(cores, 2)
    expect_false(logical)
    expect_equal(leave_free, 0)
    list(coef = model$coef)
  })
  result <- validate_cie_model(model, k = 3, custom_sky_coef = model$coef,
    method = c("BFGS", "CG"), parallel = TRUE, cores = 2,
    logical = FALSE, leave_free = 0)
  expect_equal(calls, 3)
  expect_equal(result$rmse, 0, tolerance = 1e-12)
})

test_that("two observations and single-row training tables are supported", {
  model <- cie_validation_fixture()
  model$relative_radiance$sky_points <- model$relative_radiance$sky_points[1:2, ]
  result <- validate_cie_model(model, k = 2, std_sky_no = 4)
  expect_equal(result$pred, result$obs, tolerance = 1e-12)
  expect_equal(vapply(result$fold_models,
    function(x) nrow(x$relative_radiance$sky_points), integer(1)), c(`1` = 1L, `2` = 1L))
})
