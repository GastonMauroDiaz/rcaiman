test_that("radial geometry retains a rectangular grid and an off-center axis", {
  g <- set_internal_geometry(c(8, 6), list(
    center_of_vision = c(3, 2),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = 10,
    normalization_angle = 90
  ))
  expect_s4_class(g, "SpatRaster")
  expect_equal(dim(g), c(6, 8, 2))
  expect_identical(names(g), c("off_axis", "internal_azimuth"))
  expect_equal(unname(as.vector(terra::ext(g))), c(0, 8, 0, 6))
  expect_equal(terra::res(g), c(1, 1))
  expect_null(attr(g, "projection_coef"))
  expect_null(attr(g, "parameters"))
  expect_null(attr(g, "log"))

  xy <- terra::xyFromCell(g, seq_len(terra::ncell(g)))
  dx <- xy[, 1] - 3
  dy <- 6 - xy[, 2] - 2
  expect_equal(terra::values(g)[, 1], sqrt(dx^2 + dy^2) * 9,
               tolerance = 1e-8)
  # Four diagonal rays around the snapped center, clockwise from the top.
  cells <- terra::cellFromRowCol(g, c(2, 3, 3, 2), c(4, 4, 3, 3))
  expect_equal(terra::values(g)[cells, 2], c(45, 135, 225, 315))
})

test_that("normalization does not impose coverage or an image boundary", {
  radius <- sqrt(0.5)
  g <- set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "radial_projection",
    parameters = 1 / (pi / 3),
    normalization_radius = radius,
    normalization_angle = 60,
    max_off_axis_angle = 150
  ))
  values <- terra::values(g)
  expect_equal(unname(values[6, 1]), 60, tolerance = 1e-8)
  expect_equal(unname(values[2, 1]), 60 * sqrt(5), tolerance = 1e-8)
  expect_true(all(is.na(values[1, ])))
  expect_identical(is.na(values[, 1]), is.na(values[, 2]))

  # Same physical mapping, expressed using a different normalization.
  equivalent <- set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = radius * 1.5,
    normalization_angle = 90,
    max_off_axis_angle = 150
  ))
  expect_equal(terra::values(equivalent), values, tolerance = 1e-8)

  default_limit <- set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "radial_projection",
    parameters = 1 / (pi / 3),
    normalization_radius = radius,
    normalization_angle = 60
  ))
  expect_true(all(is.na(terra::values(default_limit)[2, ])))
})

test_that("radial inversion handles a nonlinear increasing projection", {
  angle <- pi / 3
  g <- set_internal_geometry(c(6, 4), list(
    center_of_vision = c(3, 2),
    model = "radial_projection",
    parameters = c(0, 1 / angle^2),
    normalization_radius = 2,
    normalization_angle = 60
  ))
  xy <- terra::xyFromCell(g, seq_len(terra::ncell(g)))
  relative_radius <- sqrt((xy[, 1] - 3)^2 + (4 - xy[, 2] - 2)^2) / 2
  expect_equal(terra::values(g)[, 1], 60 * sqrt(relative_radius),
               tolerance = 1e-8)
})

test_that("radial input checks do not repair coefficients or fold rays", {
  build <- function(parameters, normalization_angle = 90,
                    max_off_axis_angle = 90) {
    set_internal_geometry(c(4, 4), list(
      center_of_vision = c(2, 2),
      model = "radial_projection",
      parameters = parameters,
      normalization_radius = 2,
      normalization_angle = normalization_angle,
      max_off_axis_angle = max_off_axis_angle
    ))
  }
  expect_error(build(1), "must equal one")
  expect_error(build(numeric()), "empty")
  expect_error(build(NA_real_), "NA")
  expect_error(build(2 / pi, max_off_axis_angle = 181), "180")
  expect_error(set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "radial_projection",
    parameters = 2 / pi
  )),
               "normalization_radius")
  # Normalized at 60 degrees, but folds before the requested limit.
  normalization <- 3 * pi / 3 - 2 * (pi / 3)^2
  expect_error(build(c(3, -2) / normalization, 60, 120),
               "strictly increasing")
  expect_error(set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2.5, 2),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = 2,
    normalization_angle = 90
  )),
               "integer-like")
})

test_that("anisotropic model predicts internal clock angles on the declared grid", {
  parameters <- list(off_axis_coef = 1,
                     sin_internal_azimuth_coef = c(0, 1),
                     cos_internal_azimuth_coef = c(-1, 0), atan_scale = 0.6)
  g <- set_internal_geometry(c(8, 6), list(
    center_of_vision = c(3, 2),
    model = "anisotropic_projection",
    parameters = parameters
  ))
  values <- terra::values(g)
  cells <- terra::cellFromRowCol(g, c(2, 3, 3, 2), c(4, 4, 3, 3))
  expect_equal(values[cells, 2], c(45, 135, 225, 315))
  expect_equal(values[cells, 1], rep(atan(sqrt(0.5) / 5 / 0.6) * 180 / pi, 4))
  expect_identical(names(g), c("off_axis", "internal_azimuth"))

  # A non-isotropic angular mapping changes the ray's azimuth.
  parameters$sin_internal_azimuth_coef <- c(0, 2)
  distorted <- set_internal_geometry(c(8, 6), list(
    center_of_vision = c(3, 2),
    model = "anisotropic_projection",
    parameters = parameters
  ))
  expect_equal(unname(terra::values(distorted)[cells[1], 2]), atan(2) * 180 / pi)
  expect_equal(terra::values(distorted)[, 1], values[, 1])

  # This model uses the image semidiagonal, not radial-model normalization.
  ignored <- set_internal_geometry(c(8, 6), list(
    center_of_vision = c(3, 2),
    model = "anisotropic_projection",
    parameters = parameters,
    normalization_radius = 400,
    normalization_angle = 30
  ))
  expect_equal(terra::values(ignored), terra::values(distorted))

  parameters$sin_internal_azimuth_coef <- 1
  expect_error(set_internal_geometry(c(8, 6), list(
    center_of_vision = c(3, 2),
    model = "anisotropic_projection",
    parameters = parameters
  )),
               "parameters$sin_internal_azimuth_coef", fixed = TRUE)
})

test_that("anisotropic targets are explicit and predictions are not hemisphere-clipped", {
  parameters <- list(off_axis_coef = 2,
                     sin_internal_azimuth_coef = c(0, 1),
                     cos_internal_azimuth_coef = c(-1, 0), atan_scale = 0.1)
  g <- set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "anisotropic_projection",
    parameters = parameters
  ))
  expect_true(any(terra::values(g)[, 1] > 90))
  parameters$sin_internal_azimuth_coef <- c(0, 0)
  parameters$cos_internal_azimuth_coef <- c(0, 0)
  undefined <- set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "anisotropic_projection",
    parameters = parameters
  ))
  expect_true(all(is.na(terra::values(undefined)[, 2])))
  expect_equal(terra::values(undefined)[, 1], terra::values(g)[, 1])
  names(parameters)[1] <- "zen_coef"
  expect_error(set_internal_geometry(c(4, 4), list(
    center_of_vision = c(2, 2),
    model = "anisotropic_projection",
    parameters = parameters
  )),
               "off_axis_coef")
})
