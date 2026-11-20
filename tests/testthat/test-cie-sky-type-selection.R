sky_type_fixture <- function() {
  sun <- data.frame(zenith = pi / 4, azimuth = .3)
  samples <- expand.grid(zenith = seq(.1, 1.4, length.out = 4),
                         azimuth = seq(0, 5.5, length.out = 5))
  coef <- as.numeric(cie_table[4, 1:5])
  samples$relative_radiance <- .cie_sky_model(samples$azimuth, samples$zenith,
    sun$azimuth, sun$zenith, coef[1], coef[2], coef[3], coef[4], coef[5])
  list(rr = list(zenith_dn = 1000, sky_points = samples), sun = sun)
}

test_that("sky-type selection tries the union once in standard table order", {
  f <- sky_type_fixture()
  starts <- list()
  local_mocked_bindings(optim = function(par, fn, method, ...) {
    starts[[length(starts) + 1L]] <<- par
    list(par = par, value = fn(par), convergence = 0L)
  }, .package = "stats")
  selections <- list("Overcast", c("Clear", "Partly cloudy"),
                     c("Clear", "Clear"), NULL)
  for (selection in selections) {
    starts <- list()
    fit <- fit_cie_model(f$rr, f$sun, general_sky_type = selection, method = "BFGS")
    rows <- if (is.null(selection)) seq_len(nrow(cie_table)) else
      which(cie_table$general_sky_type %in% selection)
    expect_equal(unname(do.call(rbind, starts)), unname(as.matrix(cie_table[rows, 1:5])))
    expect_equal(nrow(fit$attempts), 2 * length(rows))
  }
  starts <- list()
  fit_cie_model(f$rr, f$sun, std_sky_no = 4,
                general_sky_type = c("Clear", "Partly cloudy"), method = "BFGS")
  expect_length(starts, 1)
  expect_equal(unname(starts[[1]]), as.numeric(cie_table[4, 1:5]))
})

test_that("sky-type selection rejects empty, unknown and excessive selectors", {
  f <- sky_type_fixture()
  for (selection in list(character(), c("Overcast", "Partly cloudy", "Clear"),
                         c("Clear", NA_character_), c("Clear", "unknown"), 1)) {
    expect_error(fit_cie_model(f$rr, f$sun, general_sky_type = selection), "general_sky_type")
  }
  expect_error(fit_cie_model(f$rr, f$sun,
    custom_sky_coef = as.numeric(cie_table[4, 1:5]),
    general_sky_type = c("Overcast", "Partly cloudy")), "combine")
})

test_that("solar refinement and validation accept and transmit two sky types", {
  f <- sky_type_fixture()
  local_mocked_bindings(optim = function(par, fn, method, ...) {
    list(par = par, value = fn(par), convergence = 0L)
  }, .package = "stats")
  model <- fit_cie_model(f$rr, f$sun, std_sky_no = 4, method = "BFGS")
  types <- c("Overcast", "Partly cloudy")
  validation <- .cie_validation_predictions(model, k = 2, fold_id = NULL,
    custom_sky_coef = NULL, std_sky_no = NULL, general_sky_type = types,
    method = "BFGS", parallel = FALSE, cores = NULL, logical = TRUE, leave_free = 1)
  n_skies <- sum(cie_table$general_sky_type %in% types)
  expect_length(validation$fold_models, 2)
  for (fit in validation$fold_models) expect_equal(nrow(fit$attempts), 2 * n_skies)

  calls <- list()
  local_mocked_bindings(fit_cie_model = function(relative_radiance, sun_angles, ...) {
    calls[[length(calls) + 1L]] <<- list(...)
    model
  })
  optimize_sun_angles(model, max_iter = 0, method_sun = NULL,
                     general_sky_type = types, method = "BFGS")
  expect_length(calls, 1)
  expect_identical(calls[[1]]$general_sky_type, types)
  expect_error(optimize_sun_angles(model, general_sky_type = character()), "one or two")
})
