# These reference tests protect against incidental algorithmic changes within
# the intended use of the package during the spatial-convention migration.

test_that("fibonacci_directions() preserves the reference angular lattice", {
  observed <- fibonacci_directions(60)

  expected <- data.frame(
    off_axis = c(
      0.4111378623223475,
      0.7227342478134157,
      0.9479697413828937,
      1.1410208954903689,
      1.3181160716528177,
      1.4873662401842815
    ),
    internal_azimuth = c(
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


test_that("fibonacci_directions() supports spherical caps beyond 90 degrees off axis", {
  observed <- fibonacci_directions(60, max_off_axis_angle = 120)

  expect_equal(nrow(observed), 9)
  expect_gt(max(observed$off_axis), pi / 2)
  expect_lte(max(observed$off_axis), 2 * pi / 3)
})


test_that("internal Fibonacci directions have a deterministic angular contract", {
  directions <- fibonacci_directions(15, max_off_axis_angle = 60)
  expect_identical(directions, fibonacci_directions(15, max_off_axis_angle = 60))
  expect_true(all(directions$off_axis > 0 & directions$off_axis < pi / 3))
  expect_true(all(directions$internal_azimuth >= 0 & directions$internal_azimuth < 2 * pi))
  expect_identical(names(fibonacci_directions(180, 1)), names(directions))
  expect_equal(nrow(fibonacci_directions(180, 1)), 0)
  expect_error(fibonacci_directions(0), "positive")
  expect_error(fibonacci_directions(Inf), "Inf")
  expect_error(fibonacci_directions(181), "180")
  expect_error(fibonacci_directions(10, max_off_axis_angle = 181), "180")
  expect_error(fibonacci_directions(10, max_off_axis_angle = NA_real_), "NA")
  neighbors <- nn_spherical(sampling_directions = directions,
    query_directions = directions, k = 1,
    angle_columns = c("off_axis", "internal_azimuth"), parallel = FALSE)
  expect_identical(vapply(neighbors, function(x) x$index, integer(1)), seq_len(nrow(directions)))
  expect_equal(vapply(neighbors, function(x) x$distance, numeric(1)), rep(0, nrow(directions)))
})

test_that("reproject_points() returns continuous target coordinates", {
  sampling_points <- data.frame(
    col = 1:5,
    row = 1:5,
    zenith = c(0, rep(pi / 4, 4)),
    azimuth = c(NA, 0, pi / 2, pi, 3 * pi / 2),
    id = letters[1:5],
    x = 99,
    y_down = 99,
    y = 99,
    relative_radius = 0.5,
    image_angle = 0
  )

  observed <- reproject_points(
    sampling_points,
    list(model = "radial_projection", center_of_vision = c(50, 50),
         parameters = equidistant_coef, normalization_radius = 50,
         normalization_angle = 90)
  )

  expect_identical(names(observed), c("x", "y_down", "zenith", "azimuth", "id"))
  expect_equal(observed$x, c(50, 50, 25, 50, 75), tolerance = 1e-12)
  expect_equal(
    observed$y_down,
    c(50, 25, 50, 75, 50),
    tolerance = 1e-12
  )
})


test_that("reproject_points() uses the declared normalization and center", {
  observed <- reproject_points(
    data.frame(zenith = pi / 3, azimuth = 3 * pi / 2),
    list(model = "radial_projection", center_of_vision = c(10, 20),
         parameters = 1 / (pi / 3), normalization_radius = 12.5,
         normalization_angle = 60)
  )

  expect_equal(observed$x, 22.5, tolerance = 1e-12)
  expect_equal(observed$y_down, 20, tolerance = 1e-12)
})
