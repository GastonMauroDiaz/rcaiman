flat_field_angles <- function(values) {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 2, xmax = 5,
                    ymin = 1, ymax = 3, crs = "")
  terra::values(r) <- values
  names(r) <- "off_axis"
  r
}

test_that("sampled directions and image predictions use the same response", {
  geometry <- list(model = "radial_projection", center_of_vision = c(4, 3),
                   parameters = 2 / pi, normalization_radius = 4,
                   normalization_angle = 90)
  g <- set_internal_geometry(c(10, 6), geometry)
  indices <- rem_unselected(data.frame(col = c(4, 5, 2, 4),
                                       row = c(3, 4, 3, 3)),
                                    !is.na(g$off_axis) & !is.na(g$internal_azimuth))
  points <- add_directions(indices, g = g)
  cases <- list(simple_polynomial = list(parameters = c(-.1, -.02)),
                free_form_polynomial = list(parameters = c(.98, -.1, -.02)),
                lang2010 = list(aperture = 5))
  for (model in names(cases)) {
    args <- c(list(model = model), cases[[model]])
    response <- do.call(calc_relative_response,
      c(list(off_axis = points$off_axis, units = "radian"), args))
    ff <- do.call(flat_field_image, c(list(off_axis = g$off_axis), args))
    cells <- terra::cellFromRowCol(ff, indices$row, indices$col)
    expect_equal(response, terra::values(ff, mat = FALSE)[cells])
    expect_equal(do.call(calc_relative_response,
      c(list(off_axis = points$off_axis * 180 / pi), args)), response)
    corrected <- apply_flat_field(100 * ff, ff)
    expect_equal(terra::values(corrected, mat = FALSE)[cells], rep(100, nrow(indices)))
  }
})

test_that("numeric evaluation retains names and permits inspection of fitted curves", {
  angle <- c(axis = 0, one = 1, missing = NA_real_)
  expect_equal(calc_relative_response(angle, "simple_polynomial", -2,
                                      units = "radian"),
               c(axis = 1, one = -1, missing = NA_real_))
  expect_equal(calc_relative_response(angle, "free_form_polynomial", 0,
                                      units = "radian"),
               c(axis = 0, one = 0, missing = NA_real_))
  expect_equal(calc_relative_response(1, "simple_polynomial",
                                      c(0, 0, 0, 0, 0, -.1), "radian"), .9)
  expect_equal(calc_relative_response(1, "free_form_polynomial",
                                      c(.98, 0, 0, 0, 0, 0, -.1), "radian"), .88)
  expect_equal(calc_relative_response(numeric(), "simple_polynomial", 0), numeric())
  expect_error(calc_relative_response(0, "ff_simple_polynomial", 0), "model")
  expect_error(calc_relative_response(0, "simple_polynomial", c(0, NA)), "NA")
  expect_error(calc_relative_response(0, "simple_polynomial", Inf), "Inf")
  expect_error(calc_relative_response(0, "lang2010", 1, aperture = 5), "NULL")
  expect_error(calc_relative_response(0, "simple_polynomial", 0, aperture = 5), "NULL")
  expect_error(calc_relative_response(pi, "lang2010", units = "radian", aperture = 5),
               "off_axis")
  expect_error(calc_relative_response(0, "simple_polynomial", 0, "degrees"), "units")
})

test_that("polynomials keep their constraints and convert degrees to radians", {
  off_axis <- flat_field_angles(c(0, 30, 60, 90, 120, NA))
  before <- terra::values(off_axis)
  simple <- flat_field_image(off_axis, "simple_polynomial", -0.1)
  free <- flat_field_image(off_axis, "free_form_polynomial", c(0.98, -0.1))
  expect_equal(terra::values(simple, mat = FALSE),
               c(1, 1 - pi / 60, 1 - pi / 30, 1 - pi / 20,
                 1 - pi / 15, NA))
  expect_equal(terra::values(simple - free, mat = FALSE),
               c(rep(0.02, 5), NA))
  expect_identical(terra::values(off_axis), before)
  expect_true(terra::compareGeom(simple, off_axis))
  expect_identical(names(simple), "flat_field")
  expect_identical(names(free), "flat_field")

})

test_that("Lang evaluates the empirical expression without renormalizing", {
  off_axis <- flat_field_angles(c(0, 30, 60, 90, NA, 0))
  ff <- flat_field_image(off_axis, "lang2010", aperture = 5)
  # Reference evaluations of the published expression, including the axis.
  expected <- c(1.020535674843593, 0.969330956796954,
                0.903201273915440, 0.688463871669368, NA, 1.020535674843593)
  expect_equal(terra::values(ff, mat = FALSE), expected, tolerance = 1e-10)
  expect_true(terra::compareGeom(ff, off_axis))
  expect_identical(names(ff), "flat_field")
})

test_that("flat fields stay with internal geometry through orientation and crop", {
  geometry <- list(model = "radial_projection", center_of_vision = c(4, 3),
                   parameters = 2 / pi, normalization_radius = 4,
                   normalization_angle = 90)
  g <- set_internal_geometry(c(10, 6), geometry)
  evaluate <- list(function(x) flat_field_image(x, "simple_polynomial", -0.1),
                   function(x) flat_field_image(x, "free_form_polynomial", c(0.98, -0.1)),
                   function(x) flat_field_image(x, "lang2010", aperture = 5))
  for (model in evaluate) {
    expected <- model(g[["off_axis"]])
    d <- direction_image(g, sensor_rotation = 90, elevation = 32.5)
    expect_identical(names(d), c("zenith", "azimuth"))
    expect_equal(terra::values(model(g[["off_axis"]])),
                 terra::values(expected))
    cropped <- crop_caim(g, c(2, 2), width = 6, height = 3)
    observed <- model(cropped[["off_axis"]])
    reference <- crop_caim(expected, c(2, 2), width = 6, height = 3)
    expect_true(terra::compareGeom(observed, reference))
    expect_equal(terra::values(observed), terra::values(reference))
    missing <- g[["off_axis"]] * NA_real_
    expect_true(all(is.na(terra::values(model(missing)))))
    expect_error(model(g), "single layer")
  }
})

test_that("flat-field arguments respect angular and formula domains", {
  off_axis <- flat_field_angles(rep(30, 6))
  for (invalid in c(-1, 181, Inf)) {
    bad <- off_axis
    terra::values(bad) <- invalid
    expect_error(flat_field_image(bad, "simple_polynomial", 0), "off_axis")
  }
  too_wide <- off_axis
  terra::values(too_wide) <- 91
  expect_error(flat_field_image(too_wide, "lang2010", aperture = 5), "off_axis")
  expect_error(flat_field_image(off_axis, "lang2010", aperture = 22.87), "aperture")
  expect_error(flat_field_image(off_axis, "lang2010", aperture = 0), "aperture")
  expect_error(flat_field_image(off_axis, "simple_polynomial", rep(0, 7)), "length")
  expect_error(flat_field_image(off_axis, "free_form_polynomial", rep(0, 8)), "length")
  expect_error(flat_field_image(off_axis, "simple_polynomial", -10), "strictly positive")
  expect_error(flat_field_image(off_axis, "free_form_polynomial", 0), "strictly positive")
  expect_warning(flat_field_image(off_axis, "simple_polynomial", 1), "upper bound")
})
