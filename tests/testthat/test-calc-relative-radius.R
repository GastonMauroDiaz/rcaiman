test_that("calc_relative_radius() preserves order during extrapolation", {
  angle <- c(pi, pi / 6, 2 * pi / 3, pi / 3)

  observed <- calc_relative_radius(
    angle,
    equal_area_coef,
    units = "radian"
  )
  expected <- vapply(
    angle,
    calc_relative_radius,
    numeric(1),
    projection_coef = equal_area_coef,
    units = "radian"
  )

  expect_equal(observed, expected)
})

test_that("forward evaluation retains the polynomial beyond 90 degrees", {
  angle <- c(edge = 120, axis = 0, reference = 60, missing = NA, middle = 90)
  coefficients <- normalize_projection_coef(c(0, 1), 60)
  expect_equal(calc_relative_radius(angle, coefficients), (angle / 60)^2)
  expect_equal(calc_relative_radius(angle * pi / 180, coefficients, "radian"),
               (angle / 60)^2)
  # A forward evaluation can inspect coefficients before accepting them as
  # a lens projection; it neither renormalizes nor changes their shape.
  expect_equal(calc_relative_radius(120, c(1, -1)),
               2 * pi / 3 - (2 * pi / 3)^2)
})

test_that("inverse evaluation recovers analytic angles and the interval ends", {
  coefficients <- normalize_projection_coef(c(0, 1), 60)
  rho <- c(edge = 4, outside = 4.1, axis = 0, reference = 1,
           missing = NA, intermediate = .25)
  expected <- c(edge = 120, outside = NA, axis = 0, reference = 60,
                missing = NA, intermediate = 30)
  expect_equal(calc_off_axis_angle(rho, coefficients, 120), expected,
               tolerance = 1e-9)
  expect_equal(calc_off_axis_angle(rho, coefficients, 120, "radian"),
               expected * pi / 180, tolerance = 1e-11)
  expect_identical(calc_off_axis_angle(c(0, 2), 2 / pi, 180), c(0, 180))
  expect_true(is.na(calc_off_axis_angle(4, coefficients)))
})

test_that("both directions agree for standard projections and changed references", {
  for (coefficients in list(equidistant_coef, equal_area_coef, orthographic_coef)) {
    for (reference in c(60, 90)) {
      parameters <- normalize_projection_coef(coefficients, reference)
      angles <- c(90, 89.999, 0, 57.5, 5, 75)
      rho <- calc_relative_radius(angles, parameters)
      expect_equal(calc_off_axis_angle(rho, parameters), angles, tolerance = 1e-9)
    }
  }
})

test_that("numeric inversion and raster geometry use the same projection", {
  geometry <- list(model = "radial_projection", center_of_vision = c(2, 3),
                   parameters = normalize_projection_coef(c(.7, -.03, .02), 60),
                   normalization_radius = 1.5, normalization_angle = 60,
                   max_off_axis_angle = 120)
  g <- set_internal_geometry(c(6, 4), geometry)
  xy <- terra::xyFromCell(g, seq_len(terra::ncell(g)))
  rho <- sqrt((xy[, 1] - 2)^2 + (4 - xy[, 2] - 3)^2) / 1.5
  expected <- calc_off_axis_angle(rho, geometry$parameters, 120)
  expect_equal(as.vector(terra::values(g$off_axis)), expected)
  expect_true(anyNA(expected))
  expect_true(any(expected > 90, na.rm = TRUE))
  expect_identical(unname(is.na(terra::values(g$off_axis))),
                    unname(is.na(terra::values(g$internal_azimuth))))
})

test_that("inversion rejects ambiguous or undefined radial models", {
  # Increasing at both ends, but decreasing on an interior interval:
  # P'(x) = (x - 0.4) * (x - 0.8).
  expect_error(calc_off_axis_angle(.1, c(.32, -.6, 1/3)), "strictly increasing")
  expect_error(calc_off_axis_angle(.1, -1), "strictly increasing")
  expect_error(calc_off_axis_angle(.1, 0), "strictly increasing")
  expect_error(calc_off_axis_angle(.1, .Machine$double.xmax, 180), "finite positive")
  for (limit in list(0, -1, 181, NA_real_, Inf, c(60, 90))) {
    expect_error(calc_off_axis_angle(.5, equidistant_coef, limit))
  }
})

test_that("evaluators handle missing values and reject malformed inputs", {
  expect_identical(calc_relative_radius(numeric(), equidistant_coef), numeric())
  expect_identical(calc_off_axis_angle(numeric(), equidistant_coef), numeric())
  expect_identical(calc_relative_radius(c(NA_real_, NA_real_), equidistant_coef),
                    c(NA_real_, NA_real_))
  expect_identical(calc_off_axis_angle(c(NA_real_, NA_real_), equidistant_coef),
                    c(NA_real_, NA_real_))
  for (coefficients in list(numeric(), NA_real_, Inf, "1")) {
    expect_error(calc_relative_radius(30, coefficients))
    expect_error(calc_off_axis_angle(.5, coefficients))
  }
  expect_error(calc_relative_radius(-1, equidistant_coef), "non-negative")
  expect_error(calc_relative_radius(181, equidistant_coef), "180")
  expect_error(calc_relative_radius(pi + .1, equidistant_coef, "radian"), "180")
  expect_error(calc_relative_radius(Inf, equidistant_coef), "Inf")
  expect_error(calc_off_axis_angle(-1, equidistant_coef), "non-negative")
  expect_error(calc_off_axis_angle(Inf, equidistant_coef), "Inf")
  expect_error(calc_off_axis_angle(.5, equidistant_coef, units = "turns"))
})
