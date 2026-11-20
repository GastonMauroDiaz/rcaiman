orientation_sky <- function(zenith, azimuth) {
  100 + 20 * cos(zenith) + 15 * sin(zenith)^2 * cos(2 * azimuth) +
    10 * sin(zenith)^3 * sin(3 * azimuth)
}

orientation_fixture <- function(
    initial = list(sensor_rotation = 0, azimuth_rotation = 0, elevation = 90),
    known = list(sensor_rotation = 20, azimuth_rotation = 22, elevation = 87),
    scale = 1) {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(80, 80), normalization_radius = 80,
    normalization_angle = 90)
  g <- set_internal_geometry(c(160, 160), geometry)
  d <- do.call(direction_image, c(list(g), initial))
  reference <- orientation_sky(d$zenith * pi/180, d$azimuth * pi/180)
  samples <- fibonacci_directions(8)
  samples <- samples[samples$off_axis < 55 * pi/180, ]
  external <- do.call(external_directions_from_internal, c(list(samples), known))
  samples$Blue <- scale * orientation_sky(external$zenith, external$azimuth)
  list(samples = samples, reference = reference, geometry = geometry,
       initial = initial, known = known, external = external, scale = scale)
}

fit_orientation_fixture <- function(f, ...) {
  optimize_orientation_by_radiometry(f$samples, "Blue", f$reference, f$geometry,
    initial_orientation = f$initial, reference_orientation = f$initial,
    reference_scale = f$scale, ...)
}

test_that("radiometry recovers a known rotation near the vertical", {
  f <- orientation_fixture(scale = 2.3)
  before <- f$samples
  fit <- fit_orientation_fixture(f)
  expect_lt(fit$metric, fit$initial_metric / 1000)
  error <- calc_spherical_distance(fit$samples$zenith, fit$samples$azimuth,
                                   f$external$zenith, f$external$azimuth)
  expect_lt(max(error) * 180/pi, 0.05)
  expect_lte(fit$rotation_angle, 10)
  expect_equal(fit$metric, mean(fit$samples$residual^2))
  expect_equal(fit$samples$residual, f$samples$Blue - fit$samples$predicted)
  expect_identical(f$samples, before)
  expect_identical(rownames(fit$samples), rownames(before))
  reconstructed <- do.call(external_directions_from_internal,
                           c(list(f$samples), fit$orientation))
  expect_equal(reconstructed$zenith, fit$samples$zenith, tolerance = 1e-12)
  expect_equal(reconstructed$azimuth, fit$samples$azimuth, tolerance = 1e-12)
  # Audit reference coordinates through the public conversion functions.
  internal <- do.call(internal_directions_from_external,
                      c(list(reconstructed), f$initial))
  coord <- coord_from_internal_directions(internal, f$geometry)
  expect_equal(coord$x, fit$samples$reference_x, tolerance = 1e-10)
  expect_equal(coord$y_down, fit$samples$reference_y_down, tolerance = 1e-10)
  expect_equal(nrow(fit$attempts), 7)
})

test_that("the reference orientation and radiometric scale remain fixed", {
  f <- orientation_fixture(
    initial = list(sensor_rotation = 15, azimuth_rotation = 60, elevation = 32.5),
    known = list(sensor_rotation = 18, azimuth_rotation = 62, elevation = 34),
    scale = 0.7)
  fit <- fit_orientation_fixture(f)
  error <- calc_spherical_distance(fit$samples$zenith, fit$samples$azimuth,
                                   f$external$zenith, f$external$azimuth)
  expect_lt(max(error) * 180/pi, 0.05)
  expect_lt(fit$metric, fit$initial_metric / 1000)
  # Pre-scaling the reference is equivalent to the explicit multiplier.
  f$reference <- f$reference * f$scale
  f$scale <- 1
  scaled <- fit_orientation_fixture(f)
  expect_equal(scaled$samples$predicted, fit$samples$predicted, tolerance = 1e-6)
})

test_that("constant skies preserve the initial orientation and do not fit gain", {
  f <- orientation_fixture()
  terra::values(f$reference) <- 100
  f$samples$Blue <- 250
  f$scale <- 2
  set.seed(753)
  rng <- .Random.seed
  fit <- fit_orientation_fixture(f, n_starts = 9, max_iter = 30)
  expect_identical(.Random.seed, rng)
  expect_equal(fit$metric, 2500)
  expect_equal(fit$initial_metric, 2500)
  expect_equal(fit$rotation_angle, 0)
  expect_equal(fit$orientation, f$initial)
  expect_equal(fit, fit_orientation_fixture(f, n_starts = 9, max_iter = 30))
})

test_that("bilinear lookup respects centers, boundaries, and missing support", {
  pixels <- matrix(c(1, 3, 5, 7, 9, 11), nrow = 2, byrow = TRUE)
  expect_equal(.orientation_bilinear(pixels, c(.5, 1, 1.5, 2.5),
                                    c(.5, 1, 1.5, 1.5)), c(1, 5, 9, 11))
  expect_true(all(is.na(.orientation_bilinear(pixels, c(0, 3, NA), c(1, 1, 1)))))
  pixels[2, 2] <- NA
  expect_equal(.orientation_bilinear(pixels, .5, .5), 1)
  expect_true(is.na(.orientation_bilinear(pixels, 1, 1)))
  expect_equal(.orientation_bilinear(matrix(8, 1, 1), .5, .5), 8)
  expect_identical(.orientation_bilinear(matrix(13.7, 4, 4),
    c(.6, 1.234, 2.987), c(1.111, 3.141, .823)), rep(13.7, 3))
})

test_that("support and argument failures are explicit", {
  f <- orientation_fixture()
  expect_error(fit_orientation_fixture(f, max_rotation = 181), "at most 180")
  expect_error(fit_orientation_fixture(f, max_rotation = 0), "positive")
  expect_error(fit_orientation_fixture(f, n_starts = 0), "positive")
  expect_error(fit_orientation_fixture(f, method = "SANN"), "arg")
  f$samples$Blue[2] <- NA
  expect_error(fit_orientation_fixture(f), "NA")
  f$samples$Blue[2] <- 100
  f$samples$internal_azimuth[2] <- NA
  expect_error(fit_orientation_fixture(f), "Azimuth")
  f <- orientation_fixture()
  coord <- coord_from_internal_directions(f$samples, f$geometry)
  indices <- indices_from_coord(coord$x, coord$y_down)
  cells <- terra::cellFromRowCol(f$reference, indices$row, indices$col)
  f$reference[cells[2]] <- NA
  expect_error(fit_orientation_fixture(f), "support for rows:.*2")
})

test_that("bounded search retains all samples and improves even with few iterations", {
  f <- orientation_fixture()
  fit <- fit_orientation_fixture(f, max_rotation = 1, max_iter = 20)
  expect_lte(fit$rotation_angle, 1)
  expect_lte(fit$metric, fit$initial_metric)
  expect_equal(nrow(fit$samples), nrow(f$samples))
  expect_true(all(is.finite(fit$samples$predicted)))
  expect_true(any(fit$attempts$convergence != 0, na.rm = TRUE))
})

test_that("starts outside coverage are reported without changing the sample set", {
  f <- orientation_fixture()
  terra::values(f$reference) <- 100
  f$samples <- internal_directions_from_coord(
    data.frame(x = c(80.5, 81.5, 79.5), y_down = rep(.5, 3)), f$geometry)
  f$samples$Blue <- 100
  fit <- fit_orientation_fixture(f, max_iter = 20)
  expect_true(any(!is.finite(fit$attempts$initial_metric)))
  skipped <- !is.finite(fit$attempts$initial_metric)
  expect_true(all(is.na(fit$attempts$convergence[skipped])))
  expect_equal(fit$rotation_angle, 0)
  expect_equal(nrow(fit$samples), 3)
  expect_equal(fit$metric, 0)
})
