sampling_geometry <- function() {
  set_internal_geometry(c(10, 6), list(
    model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(5, 3), normalization_radius = 3, normalization_angle = 90))
}

test_that("directions append radians in cell order and allow successive calls", {
  g <- sampling_geometry()
  d <- direction_image(g)
  indices <- data.frame(col = c(6, 5, 6), row = c(3, 3, 3),
                        id = factor(c("b", "a", "b")), note = c(NA, "keep", "again"),
                        row.names = c("first", "second", "third"))
  original <- indices
  both <- add_directions(indices, g, d)
  expect_named(both, c(names(indices), "off_axis", "internal_azimuth", "zenith", "azimuth"))
  expect_identical(both[names(indices)], indices)
  expect_equal(both$off_axis, rep(sqrt(.5) * pi / 6, 3))
  expect_equal(both$internal_azimuth, c(pi/4, 7*pi/4, pi/4))
  expect_equal(both$zenith, both$off_axis)
  expect_equal(both$azimuth, c(7*pi/4, pi/4, 7*pi/4))
  internal <- add_directions(indices, g = g)
  external <- add_directions(indices, d = d)
  expect_named(internal, c(names(indices), "off_axis", "internal_azimuth"))
  expect_named(external, c(names(indices), "zenith", "azimuth"))
  expect_identical(add_directions(internal, d = d), both)
  expect_equal(add_directions(external, g)[names(both)], both)
  expect_identical(indices, original)
})

test_that("the new functions compose on partial and cropped geometries", {
  g <- sampling_geometry()
  d <- direction_image(g, sensor_rotation = 90, elevation = 32.5)
  indices <- data.frame(col = c(1, 6, 7, 11), row = c(1, 3, 4, 1), id = letters[1:4])
  expect_error(add_directions(indices, g), "rem_unselected")
  expect_error(add_directions(indices, d = d), "rem_unselected")
  retained <- rem_unselected(indices, !is.na(g$off_axis) & !is.na(g$internal_azimuth))
  expect_identical(retained, indices[c(2,3), ])
  expected <- extract_dn(d, retained, use_window = FALSE)
  observed <- add_directions(retained, d = d)
  expect_equal(observed$zenith, expected$zenith * pi / 180)
  expect_equal(observed$azimuth, expected$azimuth * pi / 180)
  cropped <- crop_caim(g, upper_left_indices = c(4, 2), width = 4, height = 3)
  local <- data.frame(col = c(2, 3), row = c(2, 2))
  original <- data.frame(col = c(5, 6), row = c(3, 3))
  expect_equal(add_directions(local, cropped)[c("off_axis", "internal_azimuth")],
               add_directions(original, g)[c("off_axis", "internal_azimuth")])
})

test_that("external poles remain defined with missing azimuth", {
  d <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 4, nlyrs = 2))
  names(d) <- c("zenith", "azimuth")
  values(d) <- cbind(c(0, 180, 45, NA), rep(NA_real_, 4))
  indices <- data.frame(col = 1:4, row = 1)
  expect_equal(rem_unselected(indices, !is.na(d$zenith)), indices[1:3, ])
  expect_equal(nrow(rem_unselected(indices, !is.na(d$zenith) & !is.na(d$azimuth))), 0)
  observed <- add_directions(indices[1:2, ], d = d)
  expect_identical(observed$zenith, c(0, pi))
  expect_true(all(is.na(observed$azimuth)))
  expect_error(add_directions(indices[3, ], d = d), "missing")
  expect_error(add_directions(indices[4, ], d = d), "missing")
})

test_that("direction sources require their named layers and matching grids", {
  g <- sampling_geometry()
  d <- direction_image(g)
  indices <- data.frame(col = 5, row = 3)
  expect_error(add_directions(indices), "Supply")
  expect_error(add_directions(indices, g[[1]]), "two layers")
  expect_error(add_directions(indices, g[[2:1]]), "in that order")
  expect_error(add_directions(indices, d = g), "zenith")
  moved <- terra::shift(d, dx = 1)
  expect_error(add_directions(indices, g, moved), "same geometry")
  bad <- terra::deepcopy(g)
  angles <- values(bad)
  angles[25, 1] <- 181
  values(bad) <- angles
  expect_error(add_directions(indices, bad), "degrees")
  bad <- terra::deepcopy(g)
  angles <- values(bad)
  angles[25, 2] <- NA_real_
  values(bad) <- angles
  expect_error(add_directions(indices, bad), "missing")
  # Columns belonging to the other source can be retained, but not replaced.
  internal <- add_directions(indices, g)
  expect_error(add_directions(internal, g), "reserved columns")
  external <- add_directions(indices, d = d)
  expect_error(add_directions(external, d = d), "reserved columns")
  expect_error(add_directions(data.frame(x = .5, y_down = .5), g), "col, row")
})

test_that("mask selection preserves order, repeated cells and auxiliary values", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 3, nlyrs = 2))
  names(r) <- c("col", "note")
  values(r) <- cbind(c(1, NA, 3, 4, 5, 6), c(1, 2, NA, 4, 5, 6))
  indices <- data.frame(col = c(1, 2, 3, 3, 1, 4, 1), row = c(1, 1, 1, 2, 1, 1, 3),
                        note = c(NA, "b", "c", NA, "e", "f", "g"))
  original <- indices
  m <- !is.na(r[[1]]) & !is.na(r[[2]])
  expect_identical(rem_unselected(indices, m), indices[c(1,4,5), ])
  expect_identical(rem_unselected(indices, !is.na(r[[1]])), indices[c(1,3,4,5), ])
  # Explicit exclusions differ from missing image data.
  m[] <- c(TRUE, FALSE, NA, FALSE, TRUE, FALSE)
  expect_identical(rem_unselected(indices, m), indices[c(1,5), ])
  expect_identical(names(r), c("col", "note"))
  expect_identical(indices, original)
})

test_that("empty selections retain the expected table structure", {
  g <- sampling_geometry()
  indices <- data.frame(col = integer(), row = integer(), id = character())
  expect_identical(rem_unselected(indices, !is.na(g$off_axis) & !is.na(g$internal_azimuth)), indices)
  result <- add_directions(indices, g, direction_image(g))
  expect_equal(nrow(result), 0)
  expect_named(result, c(names(indices), "off_axis", "internal_azimuth", "zenith", "azimuth"))
  expect_type(result$off_axis, "double")
  r <- terra::rast(nrows = 2, ncols = 3)
  one <- data.frame(col = 1, row = 1, id = "a")
  expect_error(rem_unselected(one, r), "logical mask")
  values(r) <- rep(FALSE, ncell(r))
  expect_identical(rem_unselected(one, r), one[FALSE, ])
  values(r) <- rep(NA, ncell(r))
  expect_identical(rem_unselected(one, r), one[FALSE, ])
})

test_that("selection requires a single logical mask", {
  r <- terra::rast(nrows = 2, ncols = 3)
  r[] <- c(0, 1, 0, 1, 0, 1)
  indices <- data.frame(col = 1:3, row = 1)
  expect_error(rem_unselected(indices, r), "logical")
  m <- r == 1
  expect_error(rem_unselected(indices, c(m, m)), "single layer")
  expect_identical(rem_unselected(indices, m), indices[2, , drop = FALSE])
  expect_error(rem_unselected(indices, TRUE), "SpatRaster")
})

test_that("index filtering rejects continuous and invalid indices", {
  g <- sampling_geometry()
  for (operation in list(function(x) rem_unselected(x, !is.na(g$off_axis) & !is.na(g$internal_azimuth)),
                        function(x) add_directions(x, g))) {
    expect_error(operation(data.frame(col = 1.5, row = 1)), "integer-like")
    expect_error(operation(data.frame(col = 1, row = 0)), "positive")
    expect_error(operation(data.frame(col = 1, row = NA_real_)), "NA")
    expect_error(operation(data.frame(x = .5, y_down = .5)), "col, row")
  }
})
