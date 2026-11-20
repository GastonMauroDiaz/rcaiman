# Known rays make the orientation tests independent of either lens model.
orientation_rays <- function(off_axis, internal_azimuth) {
  g <- set_internal_geometry(c(length(off_axis), 1), list(
    center_of_vision = c(0, 0),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = 100,
    normalization_angle = 90
  ))
  terra::values(g) <- cbind(off_axis, internal_azimuth)
  g
}

test_that("upward directions preserve internal geometry and compass convention", {
  g <- set_internal_geometry(c(12, 8), list(
    center_of_vision = c(5, 3),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = 4,
    normalization_angle = 90
  ))
  before <- terra::values(g)
  d <- direction_image(g)
  expect_identical(names(d), c("zenith", "azimuth"))
  expect_equal(dim(d), c(8, 12, 2))
  expect_true(terra::compareGeom(g, d, lyrs = FALSE))
  expect_identical(terra::values(g), before)
  expect_equal(terra::values(d)[, "zenith"], before[, 1], tolerance = 1e-12)
  expect_equal(terra::values(d)[, "azimuth"], (-before[, 2]) %% 360,
               tolerance = 1e-12)

  rays <- orientation_rays(c(0, rep(30, 4)), c(NA, 0, 90, 180, 270))
  angles <- unname(terra::values(direction_image(rays))[, c("zenith", "azimuth")])
  expect_equal(angles[, 1], c(0, rep(30, 4)))
  expect_equal(angles[, 2], c(NA, 0, 270, 180, 90))
})

test_that("hinge axis elevation is independent of portrait and heading", {
  g <- orientation_rays(0, NA_real_)
  orientations <- data.frame(sensor = c(0, 90, 180, 270, 12.5),
                             heading = c(0, 180, 270, 90, 123))
  for (i in seq_len(nrow(orientations))) {
    d <- direction_image(g, orientations$sensor[i], orientations$heading[i], 32.5)
    angles <- unname(terra::values(d)[, c("zenith", "azimuth")])
    expect_equal(angles, c(57.5, (180 + orientations$heading[i]) %% 360),
                 tolerance = 1e-12)
  }
  expect_true(is.na(terra::values(g)[1, "internal_azimuth"]))
})

test_that("clockwise sensor rotation follows the operator's view", {
  g <- orientation_rays(rep(90, 4), c(0, 90, 180, 270))
  horizontal <- unname(terra::values(direction_image(g, elevation = 0)))
  expect_equal(horizontal[, 1], c(0, 90, 180, 90))
  expect_equal(horizontal[, 2], c(NA, 270, NA, 90))

  portrait <- unname(terra::values(direction_image(g, 90, elevation = 0)))
  expect_equal(portrait[, 1], c(90, 180, 90, 0))
  expect_equal(portrait[, 2], c(270, NA, 90, NA))

  upward <- direction_image(g, sensor_rotation = 90)
  expect_equal(unname(terra::values(upward)[, "azimuth"]), c(270, 180, 90, 0))
  turned <- direction_image(g, azimuth_rotation = 90)
  expect_equal(unname(terra::values(turned)[, "azimuth"]), c(90, 0, 270, 180))
})

test_that("singular and missing rays are handled without hemisphere clipping", {
  g <- orientation_rays(c(57.5, 122.5, 60, NA, 30, 0, 180),
                         c(0, 180, 180, 90, NA, NA, NA))
  before <- terra::values(g)
  d <- direction_image(g, elevation = 32.5)
  observed <- unname(terra::values(d)[, c("zenith", "azimuth")])
  expect_equal(observed[, 1], c(0, 180, 117.5, NA, NA, 57.5, 122.5),
               tolerance = 1e-12)
  expect_equal(observed[, 2], c(NA, NA, 180, NA, NA, 180, 0),
               tolerance = 1e-12)
  expect_identical(terra::values(g), before)

  downward <- direction_image(orientation_rays(0, NA_real_), elevation = -90)
  expect_equal(unname(terra::values(downward)[, c("zenith", "azimuth")]), c(180, NA))

  missing <- orientation_rays(c(NA_real_, NA_real_), c(NA_real_, NA_real_))
  expect_true(all(is.na(terra::values(direction_image(missing)))))
})

test_that("orientation is a rigid rotation for arbitrary ray pairs", {
  g <- orientation_rays(c(12, 46, 77, 115), c(4, 127, 213, 310))
  d <- direction_image(g, sensor_rotation = 37,
                       azimuth_rotation = 123, elevation = 32.5)
  to_xyz <- function(polar, azimuth) {
    cbind(sin(polar) * cos(azimuth), sin(polar) * sin(azimuth), cos(polar))
  }
  internal_radians <- terra::values(g) * pi / 180
  external_radians <- terra::values(d) * pi / 180
  internal <- to_xyz(internal_radians[, 1], internal_radians[, 2])
  external <- to_xyz(external_radians[, 1], external_radians[, 2])
  expect_equal(tcrossprod(external), tcrossprod(internal), tolerance = 1e-12)

  wrapped <- direction_image(g, sensor_rotation = 37 - 360,
                             azimuth_rotation = 123 + 720, elevation = 32.5)
  expect_equal(terra::values(wrapped), terra::values(d), tolerance = 1e-12)
})

test_that("orientation preserves cropped grids and accepts the anisotropic model", {
  g <- set_internal_geometry(c(12, 8), list(
    center_of_vision = c(5, 3),
    model = "anisotropic_projection",
    parameters = list(off_axis_coef = 1, sin_internal_azimuth_coef = c(0, 2),
         cos_internal_azimuth_coef = c(-1, 0), atan_scale = 0.6)
  ))
  cropped <- terra::crop(g, terra::ext(2, 9, 1, 6))
  before <- terra::values(cropped)
  d <- direction_image(cropped, sensor_rotation = 270, elevation = 32.5)
  expect_true(terra::compareGeom(cropped, d, lyrs = FALSE))
  expect_identical(terra::values(cropped), before)
  full <- terra::crop(direction_image(g, 270, elevation = 32.5), cropped)
  expect_equal(terra::values(d), terra::values(full))
})

test_that("orientation validates the angular interface", {
  g <- orientation_rays(30, 45)
  expect_error(direction_image(g, elevation = 91), "elevation")
  expect_error(direction_image(g, sensor_rotation = Inf), "Inf")
  expect_error(direction_image(g, azimuth_rotation = NA_real_), "NA")
  expect_error(direction_image(g[[2:1]]), "two layers")
  expect_error(direction_image(direction_image(g)), "two layers")
  expect_error(direction_image(orientation_rays(181, 0)), "Internal angles")
  expect_error(direction_image(orientation_rays(30, 360)), "Internal angles")
})
