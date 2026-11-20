anisotropic_coord_geometry <- function() {
  list(model = "anisotropic_projection", center_of_vision = c(3, 2),
    parameters = list(off_axis_coef = c(1.1, .1),
      sin_internal_azimuth_coef = c(.1, 1.3),
      cos_internal_azimuth_coef = c(-.9, .2), atan_scale = .6))
}

test_that("anisotropic continuous predictions equal the raster at cell centers", {
  geometry <- anisotropic_coord_geometry()
  g <- set_internal_geometry(c(8, 6), geometry)
  indices <- expand.grid(col = 1:8, row = 1:6)
  coord <- coord_from_indices(indices$col, indices$row)
  calculated <- internal_directions_from_coord(coord, geometry, dim = c(8, 6))
  extracted <- add_directions(indices, g = g)
  expect_equal(calculated[c("off_axis", "internal_azimuth")],
    extracted[c("off_axis", "internal_azimuth")], tolerance = 1e-12)
  class(geometry) <- "geometry_spec"
  geometry$notes <- "Additional fields are accepted."
  expect_identical(internal_directions_from_coord(coord, geometry, c(8, 6)), calculated)
})

test_that("subpixel evaluation preserves rows and uses the declared frame scale", {
  geometry <- anisotropic_coord_geometry()
  coord <- data.frame(x = c(2.1, 2.3, 9.2, 2.1), y_down = c(1.2, 1.4, -2, 1.2),
    dn = 1:4, relative_radius = -999, image_angle = -999,
    zenith = 99, azimuth = 98, off_axis = -1, internal_azimuth = -1,
    row.names = letters[4:1])
  before <- coord
  out <- internal_directions_from_coord(coord, geometry, c(8, 6))
  retained <- setdiff(names(coord), c("off_axis", "internal_azimuth"))
  expect_identical(out[retained], coord[retained])
  expect_identical(coord, before)
  expect_identical(out$off_axis[1], out$off_axis[4])
  expect_false(isTRUE(all.equal(out$off_axis[1], out$off_axis[2])))
  expect_true(all(is.finite(out$off_axis)))
  scaled <- coord
  scaled$x <- 2 * scaled$x
  scaled$y_down <- 2 * scaled$y_down
  geometry$center_of_vision <- 2 * geometry$center_of_vision
  bigger <- internal_directions_from_coord(scaled, geometry, c(16, 12))
  expect_equal(bigger[c("off_axis", "internal_azimuth")],
    out[c("off_axis", "internal_azimuth")])
})

test_that("axis, missing coordinates, zero angular predictors and empty tables are explicit", {
  geometry <- anisotropic_coord_geometry()
  coord <- data.frame(x = c(3, NA, 1.2, 5), y_down = c(2, 1, NA, 4))
  out <- internal_directions_from_coord(coord, geometry, c(8, 6))
  expect_equal(out$off_axis[1], 0)
  expect_true(all(is.na(out$internal_azimuth[1:3])))
  expect_true(all(is.na(out$off_axis[2:3])))
  geometry$parameters$sin_internal_azimuth_coef <- c(0, 0)
  geometry$parameters$cos_internal_azimuth_coef <- c(0, 0)
  undefined <- internal_directions_from_coord(coord, geometry, c(8, 6))
  expect_true(all(is.na(undefined$internal_azimuth)))
  expect_identical(undefined$off_axis, out$off_axis)
  empty <- internal_directions_from_coord(coord[FALSE, ], geometry, c(8, 6))
  expect_identical(empty$off_axis, numeric())
  expect_identical(empty$internal_azimuth, numeric())
  expect_identical(empty[names(coord)], coord[FALSE, ])
})

test_that("anisotropic evaluation retains raw predictions and validates its parameters", {
  geometry <- anisotropic_coord_geometry()
  coord <- data.frame(x = c(1.1, 7.2), y_down = c(1.4, 3.6))
  expect_error(internal_directions_from_coord(coord, geometry), "dim")
  for (dim in list(c(8, 0), c(8, 6.2), c(8, NA), 8)) {
    expect_error(internal_directions_from_coord(coord, geometry, dim))
  }
  expect_error(internal_directions_from_coord(data.frame(x = Inf, y_down = 1), geometry, c(8, 6)), "Inf")
  geometry$parameters$off_axis_coef <- -1
  expect_true(all(internal_directions_from_coord(coord, geometry, c(8, 6))$off_axis < 0))
  geometry$parameters$off_axis_coef <- 3
  geometry$max_off_axis_angle <- 1
  expect_true(any(internal_directions_from_coord(coord, geometry, c(8, 6))$off_axis > pi / 2))
  geometry$parameters$atan_scale <- 0
  expect_error(internal_directions_from_coord(coord, geometry, c(8, 6)), "positive")
  geometry$parameters$atan_scale <- .6
  geometry$parameters$off_axis_coef <- numeric()
  expect_error(internal_directions_from_coord(coord, geometry, c(8, 6)), "empty")
  geometry$parameters$off_axis_coef <- c(1, Inf)
  expect_error(internal_directions_from_coord(coord, geometry, c(8, 6)), "Inf")
})

test_that("retired model API and its old parameter validator are absent", {
  for (name in c("add_directions_from_model", ".assert_anisotropic_projection_parameters")) {
    expect_false(exists(name, envir = asNamespace("rcaiman"), inherits = FALSE))
    expect_false(name %in% getNamespaceExports("rcaiman"))
  }
})

test_that("shared preprocessing reproduces fitted model predictions", {
  # Synthetic observations test fit/evaluation consistency, not calibration.
  calibration_coord <- data.frame(
    x = c(1, 2, 3, 5, 7, 6, 8, 2, 1, 5, 7, 4),
    y_down = c(1, 5, 2, 4, 1, 5, 3, 2, 4, 1, 6, 5),
    zen = c(1.1, 0.8, 0.2, 0.4, 0.9, 0.7, 1.2, 0.5, 1, 0.6, 1.3, 0.3),
    azi = c(0.1, 1.2, 2.3, 3.4, 4.5, 5.6, 0.8, 1.9, 3, 4.1, 5.2, 6.1)
  )
  calibration_coord <- add_image_polar(calibration_coord, 8, 6)
  atan_scale <- 0.45
  calibration_coord$atan_radius <-
    atan(calibration_coord$relative_radius / atan_scale)
  zen_fit <- stats::lm(
    zen ~ poly(atan_radius, 3, raw = TRUE) - 1,
    data = calibration_coord
  )
  sin_fit <- stats::lm(
    sin(azi) ~ sin(image_angle) + cos(image_angle) - 1,
    data = calibration_coord
  )
  cos_fit <- stats::lm(
    cos(azi) ~ sin(image_angle) + cos(image_angle) - 1,
    data = calibration_coord
  )

  prediction_coord <- add_image_polar(
    data.frame(x = c(1.5, 4.5, 6.2), y_down = c(2.8, 1.7, 4.1)),
    8, 6
  )
  parameters <- list(
    off_axis_coef = stats::coef(zen_fit),
    sin_internal_azimuth_coef = stats::coef(sin_fit),
    cos_internal_azimuth_coef = stats::coef(cos_fit),
    atan_scale = atan_scale
  )
  observed <- internal_directions_from_coord(
    prediction_coord,
    list(model = "anisotropic_projection", center_of_vision = c(4, 3),
         parameters = parameters), dim = c(8, 6)
  )
  prediction_coord$atan_radius <-
    atan(prediction_coord$relative_radius / atan_scale)

  expect_equal(
    observed$off_axis,
    unname(stats::predict(zen_fit, prediction_coord))
  )
  expect_equal(
    observed$internal_azimuth,
    unname(atan2(
      stats::predict(sin_fit, prediction_coord),
      stats::predict(cos_fit, prediction_coord)
    ) %% (2 * pi))
  )
})
