crop_test_raster <- function() {
  r <- terra::rast(nrows = 4, ncols = 5, nlyrs = 2)
  r <- set_rcaiman_geometry(r)
  terra::values(r) <- cbind(1:20, 101:120)
  names(r) <- c("first", "second")
  r
}


test_that("crop_caim() selects the requested cells exactly", {
  r <- crop_test_raster()

  observed <- crop_caim(
    r,
    upper_left_indices = c(2, 2),
    width = 3,
    height = 2
  )

  expect_equal(dim(observed), c(2, 3, 2))
  expect_identical(names(observed), names(r))
  expect_equal(
    unname(terra::values(observed)),
    cbind(c(7, 8, 9, 12, 13, 14), c(107, 108, 109, 112, 113, 114))
  )
  expect_equal(unname(as.vector(terra::ext(observed))), c(0, 3, 0, 2))
})


test_that("crop_caim() handles the full image and lower-right border", {
  r <- crop_test_raster()

  full <- crop_caim(
    r,
    upper_left_indices = c(1, 1),
    width = 5,
    height = 4
  )
  lower_right <- crop_caim(
    r,
    upper_left_indices = c(4, 3),
    width = 2,
    height = 2
  )

  expect_equal(terra::values(full), terra::values(r))
  expect_equal(
    unname(terra::values(lower_right)[, 1]),
    c(14, 15, 19, 20)
  )
})


test_that("crop_caim() validates a complete cell-aligned selection", {
  r <- crop_test_raster()

  expect_error(crop_caim(r, c(1.5, 1), 2, 2), "integer-like")
  expect_error(crop_caim(r, c(1, 1), 2.5, 2), "integer-like")
  expect_error(
    crop_caim(r, upper_left_indices = c(1, 1)),
    "must be supplied together",
    fixed = TRUE
  )
  expect_error(
    crop_caim(r, c(4, 3), width = 3, height = 2),
    "extends outside",
    fixed = TRUE
  )
})


test_that("crop_points() translates coord from the first retained cell", {
  sampling_coord <- data.frame(
    x = c(3.5, 5.2),
    y_down = c(2.5, 4.8),
    id = c("a", "b")
  )

  observed <- crop_points(
    sampling_coord,
    upper_left_indices = c(4, 3)
  )

  expect_equal(
    observed,
    data.frame(
      x = c(0.5, 2.2),
      y_down = c(0.5, 2.8),
      id = c("a", "b")
    )
  )
})


test_that("crop_points() maps retained cell centers consistently", {
  original_indices <- data.frame(
    col = c(4, 5),
    row = c(3, 5)
  )
  original_coord <- coord_from_indices(
    original_indices$col,
    original_indices$row
  )

  observed <- crop_points(original_coord, c(4, 3))
  expected <- coord_from_indices(col = c(1, 2), row = c(1, 3))

  expect_equal(observed, expected)
})


test_that("crop_points() requires coord and integer-like crop indices", {
  expect_error(
    crop_points(data.frame(col = 1, row = 1), c(1, 1)),
    "required columns: x, y_down",
    fixed = TRUE
  )
  expect_error(
    crop_points(data.frame(x = 0.5, y_down = 0.5), c(1.5, 1)),
    "integer-like"
  )
})
