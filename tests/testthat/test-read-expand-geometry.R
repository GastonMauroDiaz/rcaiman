geometry_test_image_path <- function() {
  system.file("external/DSCN4500.JPG", package = "rcaiman")
}


test_that("read_caim() centers a snapped zenith on an even crop", {
  path <- geometry_test_image_path()
  full_caim <- read_caim(path, fill_na = "none")

  observed <- read_caim(
    path,
    snapped_zenith_coord = c(5, 5),
    horizon_radius = 3,
    fill_na = "none"
  )
  expected <- crop_caim(
    full_caim,
    upper_left_indices = c(3, 3),
    width = 6,
    height = 6
  )

  expect_equal(dim(observed), c(6, 6, 3))
  expect_equal(terra::values(observed), terra::values(expected))
})


test_that("read_caim() validates its declared processing geometry", {
  path <- geometry_test_image_path()

  expect_error(
    read_caim(path, snapped_zenith_coord = c(5, 5)),
    "must be supplied together",
    fixed = TRUE
  )
  expect_error(
    read_caim(path, horizon_radius = 3),
    "must be supplied together",
    fixed = TRUE
  )
  expect_error(
    read_caim(path, snapped_zenith_coord = c(5.5, 5), horizon_radius = 3),
    "integer-like"
  )
  expect_error(
    read_caim(path, snapped_zenith_coord = c(5, 5), horizon_radius = 2.5),
    "integer-like"
  )
})


test_that("expand_non_circular() preserves cells around the snapped zenith", {
  caim <- crop_caim(
    read_caim(geometry_test_image_path()),
    upper_left_indices = c(1000, 800),
    width = 6,
    height = 4
  )
  z <- zenith_image(8, equidistant_coef)

  observed <- expand_non_circular(
    caim,
    z,
    snapped_zenith_coord = c(3, 2)
  )
  target_cells <- terra::cellFromRowCol(
    z,
    row = rep(3:6, each = 6),
    col = rep(2:7, times = 4)
  )

  expect_equal(dim(observed), c(8, 8, 3))
  expect_equal(
    unname(terra::values(observed)[target_cells, ]),
    unname(terra::values(caim))
  )
  expect_true(all(is.na(terra::values(observed)[1, ])))
})


test_that("expand_non_circular() requires a snapped grid intersection", {
  caim <- crop_caim(
    read_caim(geometry_test_image_path()),
    upper_left_indices = c(1000, 800),
    width = 6,
    height = 4
  )
  z <- zenith_image(8, equidistant_coef)

  expect_error(
    expand_non_circular(caim, z, snapped_zenith_coord = c(3.5, 2)),
    "integer-like"
  )
})
