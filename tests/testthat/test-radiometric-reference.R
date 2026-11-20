test_that("reference weighting retains inverse-square estimation and handles the axis", {
  d <- c(.07, .02, .1)
  expect_equal(.inverse_distance_weights(d), (1 / d^2) / sum(1 / d^2))
  expect_equal(.inverse_distance_weights(d * 1e-200), .inverse_distance_weights(d))
  expect_equal(.inverse_distance_weights(c(0, .1, 0)), c(.5, 0, .5))
  expect_error(.inverse_distance_weights(numeric()), "empty")
  expect_error(.inverse_distance_weights(c(1, NA)), "NA")
  expect_error(.inverse_distance_weights(-1), "non-negative")
})

test_that("relative radiance uses the shared reference estimator", {
  sky <- data.frame(col = 1:4, row = 1, zenith = c(.1, .02, .07, .4),
                    azimuth = c(0, 1, 2, 3), dn = c(100, 110, 105, 80))
  fit <- calc_relative_radiance(sky, "dn")
  expected <- sum(sky$dn[1:3] / sky$zenith[1:3]^2) / sum(1 / sky$zenith[1:3]^2)
  expect_equal(as.numeric(fit$zenith_dn), expected)
  expect_equal(fit$sky_points$relative_radiance, sky$dn / expected)
  sky$zenith[2] <- 0
  fit <- calc_relative_radiance(sky, "dn")
  expect_equal(as.numeric(fit$zenith_dn), 110)
  expect_equal(calc_relative_radiance(sky, "dn", zenith_dn = 120)$sky_points$relative_radiance,
               sky$dn / 120)
})
