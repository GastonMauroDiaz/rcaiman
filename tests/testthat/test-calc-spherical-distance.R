test_that("spherical distance compares paired direction vectors", {
  observed <- calc_spherical_distance(
    zen_1 = c(pi / 2, pi / 2, pi / 4),
    azi_1 = c(0, 0, pi / 3),
    zen_2 = c(pi / 2, pi / 2, pi / 4),
    azi_2 = c(pi / 2, pi, pi / 3)
  )

  expect_equal(observed, c(pi / 2, pi, 0), tolerance = 1e-7)
})

test_that("spherical distance broadcasts either scalar direction", {
  zen <- c(0, pi / 2, pi / 2)
  azi <- c(0, 0, pi)

  observed_many_to_one <- calc_spherical_distance(zen, azi, 0, 0)
  observed_one_to_many <- calc_spherical_distance(0, 0, zen, azi)

  expect_equal(observed_many_to_one, c(0, pi / 2, pi / 2))
  expect_equal(observed_one_to_many, observed_many_to_one)
})

test_that("spherical distance propagates missing paired directions", {
  observed <- calc_spherical_distance(
    zen_1 = c(0, NA_real_, pi / 2),
    azi_1 = c(0, 0, NA_real_),
    zen_2 = c(0, pi / 2, pi / 2),
    azi_2 = c(0, pi, 0)
  )

  expect_equal(observed[1], 0)
  expect_true(all(is.na(observed[2:3])))
  expect_length(calc_spherical_distance(numeric(), numeric(), 0, 0), 0)
})

test_that("spherical distance rejects incompatible direction lengths", {
  expect_error(
    calc_spherical_distance(c(0, 1), 0, 0, 0),
    "`zen_1` and `azi_1` must have the same length",
    fixed = TRUE
  )
  expect_error(
    calc_spherical_distance(0, 0, c(0, 1), 0),
    "`zen_2` and `azi_2` must have the same length",
    fixed = TRUE
  )
  expect_error(
    calc_spherical_distance(c(0, 1), c(0, 1), rep(0, 3), rep(0, 3)),
    "same length, or one set must contain exactly one direction",
    fixed = TRUE
  )
})

test_that("spherical distance rejects non-finite directions", {
  expect_error(
    calc_spherical_distance(c(0, Inf), c(0, 0), 0, 0),
    "Inf"
  )
  expect_error(
    calc_spherical_distance(0, 0, c(0, Inf), c(0, 0)),
    "Inf"
  )
})
