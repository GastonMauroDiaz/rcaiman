# These reference tests protect against incidental algorithmic changes within
# the intended use of the package during the spatial-convention migration.

test_that("fibonacci_zen_azi() preserves the reference angular lattice", {
  observed <- fibonacci_zen_azi(60)

  expected <- data.frame(
    zen = c(
      0.4111378623223475,
      0.7227342478134157,
      0.9479697413828937,
      1.1410208954903689,
      1.3181160716528177,
      1.4873662401842815
    ),
    azi = c(
      0,
      3.8832220774509323,
      1.4832588477222786,
      5.3664809251732057,
      2.9665176954445571,
      0.5665544657158941
    )
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("fibonacci_zen_azi() supports spherical caps beyond the horizon", {
  observed <- fibonacci_zen_azi(60, max_zenith_angle = 120)

  expect_equal(nrow(observed), 9)
  expect_gt(max(observed$zen), pi / 2)
  expect_lte(max(observed$zen), 2 * pi / 3)
})


test_that("fibonacci_indices() returns raster cells and their directions", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  observed <- fibonacci_indices(z, a, 60)
  directions <- zen_azi_from_indices(
    z,
    a,
    observed$col,
    observed$row
  )

  expect_identical(names(observed), c("col", "row", "zen", "azi"))
  expect_equal(
    observed[, c("col", "row")],
    data.frame(
      col = c(6, 7, 2, 8, 5, 3),
      row = c(4, 7, 5, 3, 10, 2)
    )
  )
  expect_equal(observed$zen, directions$zen)
  expect_equal(observed$azi, directions$azi)
})


test_that("fibonacci_indices() retains duplicated cell assignments", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  observed <- fibonacci_indices(z, a, 20)
  occupied_cells <- unique(observed[, c("col", "row")])

  expect_equal(nrow(observed), 44)
  expect_equal(nrow(occupied_cells), 43)
})


test_that("reproject_points() returns continuous target coordinates", {
  sampling_points <- data.frame(
    col = 1:5,
    row = 1:5,
    zen = c(0, rep(pi / 4, 4)),
    azi = c(0, 0, pi / 2, pi, 3 * pi / 2),
    id = letters[1:5],
    x = 99,
    y_down = 99,
    y = 99,
    relative_radius = 0.5,
    image_angle = 0
  )

  observed <- reproject_points(
    sampling_points,
    equidistant_coef,
    diameter = 100
  )

  expect_identical(names(observed), c("x", "y_down", "zen", "azi", "id"))
  expect_equal(observed$x, c(50, 50, 25, 50, 75), tolerance = 1e-12)
  expect_equal(
    observed$y_down,
    c(50, 25, 50, 75, 50),
    tolerance = 1e-12
  )
})


test_that("reproject_points() can restore a target orientation", {
  observed <- reproject_points(
    data.frame(zen = pi / 4, azi = 7 * pi / 4),
    equidistant_coef,
    diameter = 100,
    orientation = 45
  )

  expect_equal(observed$x, 50, tolerance = 1e-12)
  expect_equal(observed$y_down, 25, tolerance = 1e-12)
})
