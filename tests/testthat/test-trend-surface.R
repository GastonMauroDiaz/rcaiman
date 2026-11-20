trend_surface_samples <- function() {
  samples <- expand.grid(x = seq(0.5, 19.5, length.out = 9),
                         y_down = seq(0.5, 11.5, length.out = 8))
  samples$dn <- 0.95 - 0.001 * (samples$x - 10)^2 -
    0.002 * (samples$y_down - 6)^2 +
    0.0003 * (samples$x - 10) * (samples$y_down - 6)
  samples
}

trend_surface_template <- function() {
  set_rcaiman_geometry(terra::rast(nrows = 12, ncols = 20))
}

test_that("a compact fit reproduces spatial predictions at sample coordinates", {
  samples <- trend_surface_samples()
  samples$dn <- samples$dn +
    0.003 * sin(seq_len(nrow(samples)))
  for (degree in c(1, 2, 6)) {
    fit <- fit_trend_surface(samples, np = degree, return = "all")
    reference <- spatial::surf.ls(degree, samples$x, samples$y_down,
                                  samples$dn)
    expect_identical(names(fit), c("parameters", "model", "r2"))
    expect_s3_class(fit$model, "trls")
    expect_identical(names(fit$parameters), c("np", "beta", "rx", "ry"))
    expect_false(is.object(fit$parameters))
    expect_equal(fit$parameters, reference[c("np", "beta", "rx", "ry")])
    expect_equal(stats::fitted(fit$model), as.vector(stats::predict(
      reference, samples$x, samples$y_down)), tolerance = 1e-12)
    expect_equal(stats::residuals(fit$model), samples$dn - stats::fitted(fit$model))
    expect_equal(fit$r2, 1 - sum(stats::residuals(fit$model)^2) /
                   sum((samples$dn - mean(samples$dn))^2))
    expect_identical(fit_trend_surface(samples, np = degree), fit$parameters)
  }
  # Row order is a shared contract, independent of polynomial degree.
  shuffled <- samples[rev(seq_len(nrow(samples))), ]
  other <- fit_trend_surface(shuffled, np = degree, return = "all")
  expect_equal(stats::fitted(other$model), rev(stats::fitted(fit$model)), tolerance = 1e-10)
  expect_equal(fit_trend_surface(samples, col_id = 3, np = 2),
               fit_trend_surface(samples, np = 2))
})

test_that("an explicit prediction mask intersects the template mask", {
  parameters <- fit_trend_surface(trend_surface_samples(), np = 2)
  r <- trend_surface_template()
  reference <- trend_surface_image(r, parameters)
  m <- r
  terra::values(m) <- rep(c(TRUE, FALSE), length.out = ncell(r))
  terra::values(r) <- 1
  r[c(1, 2)] <- NA_real_
  before <- terra::values(m)
  observed <- trend_surface_image(r, parameters, m)
  expected <- terra::values(reference, mat = FALSE)
  expected[!terra::values(m, mat = FALSE)] <- NA_real_
  expected[c(1, 2)] <- NA_real_
  expect_equal(terra::values(observed, mat = FALSE), expected)
  expect_identical(terra::values(m), before)
  expect_true(terra::compareGeom(observed, r))
  terra::values(m) <- FALSE
  expect_true(all(is.na(terra::values(trend_surface_image(r, parameters, m)))))
  terra::values(m) <- TRUE
  expect_equal(terra::values(trend_surface_image(r, parameters, m)),
               terra::values(trend_surface_image(r, parameters)))
  # A mask works with a template that has no values, too.
  empty <- trend_surface_template()
  expect_equal(terra::values(trend_surface_image(empty, parameters, m)),
               terra::values(reference))
  m[1] <- NA
  expect_error(trend_surface_image(r, parameters, m), "NA-free")
  terra::values(m) <- 1
  expect_error(trend_surface_image(r, parameters, m), "logical")
  terra::values(m) <- TRUE
  terra::ext(m) <- terra::ext(1, 21, 0, 12)
  expect_error(trend_surface_image(r, parameters, m), "geometry")
})

test_that("evaluation uses pixel centers and y_down without renormalization", {
  samples <- trend_surface_samples()
  parameters <- fit_trend_surface(samples, np = 2)
  before <- parameters
  r <- trend_surface_template()
  ff <- trend_surface_image(r, parameters)
  indices <- expand.grid(col = seq_len(ncol(r)), row = seq_len(nrow(r)))
  coord <- coord_from_indices(indices$col, indices$row)
  expected <- 0.95 - 0.001 * (coord$x - 10)^2 -
    0.002 * (coord$y_down - 6)^2 +
    0.0003 * (coord$x - 10) * (coord$y_down - 6)
  expect_equal(terra::values(ff, mat = FALSE), expected, tolerance = 1e-12)
  expect_identical(names(ff), "prediction")
  expect_true(terra::compareGeom(r, ff))
  expect_identical(parameters, before)
  expect_equal(.predict_trend_surface(parameters, 10, 6), 0.95)
  corrected <- apply_flat_field(100 * ff, ff)
  expect_equal(terra::values(corrected, mat = FALSE), rep(100, ncell(r)))
  cropped <- crop_caim(ff, c(3, 2), width = 10, height = 8)
  corrected_crop <- apply_flat_field(
    crop_caim(100 * ff, c(3, 2), width = 10, height = 8), cropped
  )
  expect_equal(terra::values(corrected_crop, mat = FALSE), rep(100, 80))
})

test_that("compact parameters retain scaling and do not depend on previous fits", {
  samples <- trend_surface_samples()
  parameters <- fit_trend_surface(samples, np = 6)
  r <- trend_surface_template()
  expected <- trend_surface_image(r, parameters)
  other <- samples
  other$x <- other$x * 10 + 100
  other$y_down <- other$y_down * 5 - 50
  fit_trend_surface(other, np = 2)
  # Rebuild as plain data, as after storing the coefficients separately.
  restored <- list(np = parameters$np, beta = as.numeric(parameters$beta),
                    rx = as.numeric(parameters$rx), ry = as.numeric(parameters$ry))
  expect_equal(terra::values(trend_surface_image(r, restored)),
               terra::values(expected), tolerance = 1e-12)
})

test_that("templates mask cells without changing the fitted surface", {
  parameters <- fit_trend_surface(trend_surface_samples(), np = 2)
  r <- trend_surface_template()
  reference <- trend_surface_image(r, parameters)
  terra::values(r) <- seq_len(ncell(r))
  r[c(1, 15, 100)] <- NA_real_
  before <- terra::values(r)
  ff <- trend_surface_image(r, parameters)
  expected <- terra::values(reference)
  expected[c(1, 15, 100), ] <- NA_real_
  expect_equal(terra::values(ff), expected)
  expect_identical(terra::values(r), before)
  expect_true(all(is.na(terra::values(trend_surface_image(r * NA, parameters)))))
})

test_that("constant responses give a usable fit with undefined r2", {
  samples <- trend_surface_samples()
  samples$dn <- 0.9
  fit <- fit_trend_surface(samples, np = 2, return = "all")
  expect_true(is.na(fit$r2))
  expect_equal(stats::fitted(fit$model), rep(0.9, nrow(samples)), tolerance = 1e-12)
  expect_equal(stats::residuals(fit$model), rep(0, nrow(samples)), tolerance = 1e-12)
})

test_that("general surfaces model negative values without flat-field rules", {
  samples <- trend_surface_samples()
  samples$dn <- samples$x - 10
  parameters <- fit_trend_surface(samples, np = 1)
  r <- trend_surface_template()
  surface <- trend_surface_image(r, parameters)
  expect_equal(terra::values(surface, mat = FALSE),
               rep(seq(0.5, 19.5) - 10, 12), tolerance = 1e-12)
  terra::values(r) <- 100
  expect_error(apply_flat_field(r, surface), "strictly positive")
  terra::values(surface) <- Inf
  expect_error(apply_flat_field(r, surface), "finite")
  terra::values(surface) <- 0
  expect_error(apply_flat_field(r, surface), "strictly positive")
})

test_that("unsupported fits fail instead of returning a corrupt calibration", {
  samples <- trend_surface_samples()
  expect_error(fit_trend_surface(samples[1:5, ], np = 2), "at least 6")
  degenerate <- samples
  degenerate$y_down <- degenerate$x
  expect_error(fit_trend_surface(degenerate, np = 2), "full-rank")
  degenerate$x <- 0
  expect_error(fit_trend_surface(degenerate, np = 2), "positive width")
  expect_error(fit_trend_surface(samples, np = 7), "np")
  expect_error(fit_trend_surface(samples, return = "model"), "return")
  expect_error(fit_trend_surface(samples, col_id = 4), "col_id")
  for (bad in c(NA_real_, Inf)) {
    invalid <- samples
    invalid$dn[1] <- bad
    expect_error(fit_trend_surface(invalid, np = 2))
  }
  invalid <- samples
  invalid$x[1] <- NA_real_
  expect_error(fit_trend_surface(invalid, np = 2), "NA")
})

test_that("evaluation validates parameters, grids and finite predictions", {
  parameters <- fit_trend_surface(trend_surface_samples(), np = 2)
  r <- trend_surface_template()
  expect_error(trend_surface_image(r, parameters[-1]), "parameters")
  invalid <- parameters
  invalid$beta <- 1
  expect_error(trend_surface_image(r, invalid), "beta")
  invalid <- parameters
  invalid$ry <- rev(invalid$ry)
  expect_error(trend_surface_image(r, invalid), "increasing")
  invalid <- parameters
  invalid$beta[] <- 0
  expect_equal(terra::values(trend_surface_image(r, invalid), mat = FALSE),
               rep(0, ncell(r)))
  invalid$beta[1] <- Inf
  expect_error(trend_surface_image(r, invalid), "Inf")
  shifted <- r
  terra::ext(shifted) <- terra::ext(1, 21, 0, 12)
  expect_error(trend_surface_image(shifted, parameters), "image grid")
  expect_error(trend_surface_image(c(r, r), parameters), "single layer")
})
