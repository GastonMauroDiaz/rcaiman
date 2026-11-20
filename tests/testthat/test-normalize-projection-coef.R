test_that("normalization preserves the pixel geometry of standard projections", {
  standards <- list(equidistant_coef, equal_area_coef, orthographic_coef)
  for (coefficients in standards) {
    geometry <- list(model = "radial_projection", center_of_vision = c(10, 10),
                     parameters = coefficients, normalization_radius = 10,
                     normalization_angle = 90)
    original <- set_internal_geometry(c(20, 20), geometry)
    for (angle in c(30, 60, 90)) {
      scale <- calc_relative_radius(angle, coefficients)
      geometry$parameters <- normalize_projection_coef(coefficients, angle)
      geometry$normalization_angle <- angle
      geometry$normalization_radius <- 10 * scale
      normalized <- set_internal_geometry(c(20, 20), geometry)
      expect_equal(terra::values(normalized), terra::values(original),
                     tolerance = 1e-9)
      expect_equal(calc_relative_radius(angle, geometry$parameters), 1,
                     tolerance = 1e-14)
    }
  }
})

test_that("normalization evaluates the supplied polynomial above 90 degrees", {
  # P(theta) = theta^2, with a known normalization at 120 degrees.
  coefficients <- c(linear = 0, quadratic = 1)
  normalized <- normalize_projection_coef(coefficients, 120)
  expect_equal(normalized, c(linear = 0, quadratic = (3 / (2 * pi))^2))
  expect_equal(unname(.calc_relative_radius(2 * pi / 3, normalized)), 1)
  expect_identical(names(normalized), names(coefficients))
  expect_identical(coefficients, c(linear = 0, quadratic = 1))
  # A previous normalization at any other angle is immaterial.
  via_60 <- normalize_projection_coef(
    normalize_projection_coef(coefficients, 60), 120
  )
  expect_equal(via_60, normalized, tolerance = 1e-14)
})

test_that("normalization rejects an undefined reference", {
  for (coefficients in list(numeric(), NA_real_, Inf, "1")) {
    expect_error(normalize_projection_coef(coefficients, 60))
  }
  for (angle in list(0, -1, 181, NA_real_, Inf, c(30, 60))) {
    expect_error(normalize_projection_coef(equidistant_coef, angle))
  }
  expect_error(normalize_projection_coef(0, 60), "finite positive")
  expect_error(normalize_projection_coef(-1, 60), "finite positive")
  expect_error(normalize_projection_coef(.Machine$double.xmax, 180),
                 "finite positive")
})
