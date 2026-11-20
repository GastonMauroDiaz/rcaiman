cie_test_d <- function(zenith, azimuth) {
  d <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = length(zenith), nlyrs = 2))
  terra::values(d) <- cbind(zenith, azimuth)
  names(d) <- c("zenith", "azimuth")
  d
}

test_that("CIE images agree with the model for representative standard skies", {
  z <- c(0, 10, 40, 70, 90)
  a <- c(90, 0, 45, 180, 315)
  d <- cie_test_d(z, a)
  sun <- data.frame(zenith = 45 * pi / 180, azimuth = 30 * pi / 180)
  # All 15 coefficient sets are exercised by the numeric predictor tests.
  for (i in c(1, 8, 12)) {
    coef <- as.numeric(cie_table[i, 1:5])
    expected <- .cie_sky_model(a * pi / 180, z * pi / 180,
      sun[[2]], sun[[1]],
      coef[1], coef[2], coef[3], coef[4], coef[5])
    result <- cie_image(d, sun, coef)
    expect_equal(as.vector(result[]), expected)
    expect_equal(result[1][[1]], 1, tolerance = 1e-12)
  }
  expect_true(terra::compareGeom(result, d[[1]]))
  expect_identical(names(result), "Relative radiance or luminance")
})

test_that("analytic gradation and uniform skies are zenith-normalized", {
  d <- cie_test_d(c(0, 30, 60, 90), c(NA, 20, 170, 300))
  sun <- data.frame(zenith = 60 * pi / 180, azimuth = 40 * pi / 180)
  expect_equal(as.vector(cie_image(d, sun, c(0, 0, 0, 0, 0))[]), rep(1, 4))
  expected <- c((1 + 4 * exp(-.7 / cos(c(0, 30, 60) * pi / 180))) /
                 (1 + 4 * exp(-.7)), 1 / (1 + 4 * exp(-.7)))
  expect_equal(as.vector(cie_image(d, sun, c(4, -.7, 0, -1, 0))[]), expected)
})

test_that("missing angles, below-horizon directions and masks retain NA", {
  d <- cie_test_d(c(0, 45, 90, 91, 180, NA, 30), c(NA, NA, 90, 90, NA, 0, 180))
  sun <- data.frame(zenith = 40 * pi / 180, azimuth = 0 * pi / 180)
  coef <- as.numeric(cie_table[4, 1:5])
  sky <- cie_image(d, sun, coef)
  expect_equal(as.vector(is.na(sky[])), c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(sky[1][[1]], 1)
  mask <- d$zenith >= 0
  mask[] <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  masked <- cie_image(d, sun, coef, mask)
  expect_true(is.na(masked[1][[1]]))
  expect_equal(masked[3][[1]], sky[3][[1]])
  mask[] <- FALSE
  expect_true(all(is.na(cie_image(d, sun, coef, mask)[])))
  d[] <- NA_real_
  expect_true(all(is.na(cie_image(d, sun, coef)[])))
})

test_that("inclined partial images do not need to contain zenith", {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(40, 30), normalization_radius = 100, normalization_angle = 90)
  g <- set_internal_geometry(c(80, 60), geometry)
  d <- direction_image(g, elevation = 32.5)
  sun <- data.frame(zenith = 45 * pi / 180, azimuth = 0 * pi / 180)
  coef <- as.numeric(cie_table[4, 1:5])
  sky <- cie_image(d, sun, coef)
  expect_gt(min(d$zenith[], na.rm = TRUE), 0)
  expect_true(any(is.finite(sky[])))
  expect_true(terra::compareGeom(sky, g[[1]]))
  # Reusing external directions on another grid preserves their predictions.
  values <- as.vector(sky[])
  z <- as.vector(d$zenith[])
  a <- as.vector(d$azimuth[])
  rows <- c(1, 450, 2000, 3200)
  expect_equal(as.vector(cie_image(cie_test_d(z[rows], a[rows]), sun, coef)[]), values[rows])
  expect_null(attr(sky, "projection_coef"))
})

test_that("invalid model predictions become NA without adjusting coefficients", {
  d <- cie_test_d(c(0, 45, 90), c(0, 90, 180))
  sun <- data.frame(zenith = 40 * pi / 180, azimuth = 0 * pi / 180)
  sky <- cie_image(d, sun, c(-3, -1, 0, -1, 0))
  expect_equal(sky[1][[1]], 1)
  expect_true(all(is.na(sky[2:3][[1]])))
  expect_true(all(is.na(cie_image(d, sun, c(-1, 0, 0, -1, 0))[])))
  expect_true(all(is.na(cie_image(d, sun, c(1, 1000, 0, -1, 0))[])))
})

test_that("CIE inputs use the external raster and solar-angle contracts", {
  d <- cie_test_d(c(10, 60), c(0, 90))
  sun <- data.frame(zenith = 45 * pi / 180, azimuth = 0 * pi / 180)
  coef <- as.numeric(cie_table[4, 1:5])
  expect_error(cie_image(d[[1]], sun, coef), "two layers")
  expect_error(cie_image(d, unlist(sun), coef), "data frame")
  expect_error(cie_image(d, as.matrix(sun), coef), "data frame")
  expect_error(cie_image(d, sun[FALSE, ], coef), "one row")
  expect_error(cie_image(d, rbind(sun, sun), coef), "one row")
  expect_error(cie_image(d, data.frame(zenith = "0.5", azimuth = 0), coef), "numeric")
  expect_error(cie_image(d, data.frame(zenith = .5, azimuth = TRUE), coef), "numeric")
  expect_error(cie_image(d, data.frame(zen_deg = 45, azi_deg = 0), coef), "zenith.*azimuth")
  expect_equal(cie_image(d, rev(sun), coef)[], cie_image(d, sun, coef)[])
  expect_error(cie_image(d, data.frame(zenith = NA_real_, azimuth = 0), coef), "NA")
  expect_error(cie_image(d, data.frame(zenith = Inf, azimuth = 0), coef), "Inf")
  expect_error(cie_image(d, data.frame(zenith = -.1, azimuth = 0), coef), "non-negative")
  expect_error(cie_image(d, data.frame(zenith = pi / 4, azimuth = 2 * pi), coef), "radians")
  expect_warning(cie_image(d, data.frame(zenith = 97 * pi / 180, azimuth = 0), coef), "96 deg")
  expect_error(cie_image(d, sun, coef[1:4]), "five|5")
  expect_error(cie_image(d, sun, c(NA, coef[-1])), "NA")
  expect_error(cie_image(d, sun, coef, m = d[[1]]), "logical")
  m <- d$zenith > 0
  m[1] <- NA
  expect_error(cie_image(d, sun, coef, m), "NA-free")
  m[] <- TRUE
  terra::ext(m) <- terra::ext(0, 20, 0, 10)
  expect_error(cie_image(d, sun, coef, m), "same geometry")
})

test_that("a marked solar cell supplies radians directly through add_directions", {
  d <- cie_test_d(c(0, 45, 60), c(0, 30, 180))
  sun <- add_directions(data.frame(col = 2, row = 1), d = d)
  expect_equal(sun[c("zenith", "azimuth")],
    data.frame(zenith = pi / 4, azimuth = pi / 6))
  coef <- as.numeric(cie_table[4, 1:5])
  expected <- .cie_sky_model(c(0, 30, 180) * pi / 180,
    c(0, 45, 60) * pi / 180, pi / 6, pi / 4,
    coef[1], coef[2], coef[3], coef[4], coef[5])
  expect_equal(as.vector(cie_image(d, sun, coef)[]), expected)
})
