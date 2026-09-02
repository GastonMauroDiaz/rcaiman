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
