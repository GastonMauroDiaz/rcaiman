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


test_that("crop_caim() accepts an equivalent centered specification", {
  r <- crop_test_raster()

  centered <- crop_caim(
    r,
    snapped_zenith_coord = c(2, 2),
    horizon_radius = 1
  )
  explicit <- crop_caim(
    r,
    upper_left_indices = c(2, 2),
    width = 2,
    height = 2
  )

  expect_equal(dim(centered), c(2, 2, 2))
  expect_equal(terra::values(centered), terra::values(explicit))
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
  expect_error(
    crop_caim(r, snapped_zenith_coord = c(2, 2)),
    "must be supplied together",
    fixed = TRUE
  )
  expect_error(
    crop_caim(
      r,
      upper_left_indices = c(1, 1),
      width = 2,
      height = 2,
      snapped_zenith_coord = c(2, 2),
      horizon_radius = 1
    ),
    "not both",
    fixed = TRUE
  )
  expect_error(
    crop_caim(
      r,
      snapped_zenith_coord = c(2.5, 2),
      horizon_radius = 1
    ),
    "integer-like"
  )
  expect_error(
    crop_caim(
      r,
      snapped_zenith_coord = c(1, 1),
      horizon_radius = 2
    ),
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


test_that("point crop transformations remove previous-frame geometry", {
  sampling_coord <- data.frame(
    id = c("a", "b"),
    x = c(3.5, 5.2),
    y_down = c(2.5, 4.8),
    col = c(4L, 6L),
    row = c(3L, 5L),
    y = c(7.5, 5.2),
    relative_radius = c(0.2, 0.5),
    image_angle = c(0, pi / 2),
    zen = c(0.3, 0.7),
    azi = c(0.4, 1.2),
    dn = c(10, 20)
  )

  observed <- crop_points(
    sampling_coord,
    upper_left_indices = c(4, 3)
  )

  expect_identical(
    names(observed),
    c("x", "y_down", "id", "zen", "azi", "dn")
  )
  expect_equal(observed$x, c(0.5, 2.2))
  expect_equal(observed$y_down, c(0.5, 2.8))
  expect_identical(observed$id, sampling_coord$id)
  expect_identical(observed$zen, sampling_coord$zen)
  expect_identical(observed$azi, sampling_coord$azi)
  expect_identical(observed$dn, sampling_coord$dn)
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


test_that("crop_points() accepts an equivalent centered specification", {
  sampling_coord <- data.frame(
    x = c(4.5, 6.2),
    y_down = c(3.5, 5.8),
    id = c("a", "b")
  )

  centered <- crop_points(
    sampling_coord,
    snapped_zenith_coord = c(6, 5),
    horizon_radius = 2
  )
  explicit <- crop_points(
    sampling_coord,
    upper_left_indices = c(5, 4)
  )

  expect_equal(centered, explicit)
})


test_that("uncrop_points() is the exact inverse of crop_points()", {
  original_coord <- data.frame(
    x = c(3.5, 5.2, 3.5),
    y_down = c(2.5, 4.8, 2.5),
    id = c("a", "b", "a")
  )

  cropped_coord <- crop_points(
    original_coord,
    upper_left_indices = c(4, 3)
  )
  cropped_coord <- add_indices(cropped_coord)
  cropped_coord <- add_xy(cropped_coord, image_height = 10)
  cropped_coord$relative_radius <- c(0.2, 0.5, 0.2)
  cropped_coord$image_angle <- c(0, pi / 2, 0)
  observed <- uncrop_points(
    cropped_coord,
    upper_left_indices = c(4, 3)
  )

  expect_equal(observed, original_coord)
})


test_that("uncrop_points() accepts a centered crop specification", {
  original_coord <- data.frame(
    x = c(4.5, 6.2),
    y_down = c(3.5, 5.8),
    id = c("a", "b")
  )

  cropped_coord <- crop_points(
    original_coord,
    snapped_zenith_coord = c(6, 5),
    horizon_radius = 2
  )
  observed <- uncrop_points(
    cropped_coord,
    snapped_zenith_coord = c(6, 5),
    horizon_radius = 2
  )

  expect_equal(observed, original_coord)
})


test_that("point crop transformations preserve a valid empty data set", {
  original_coord <- data.frame(
    x = numeric(),
    y_down = numeric(),
    id = character()
  )

  cropped_coord <- crop_points(
    original_coord,
    upper_left_indices = c(4, 3)
  )
  observed <- uncrop_points(
    cropped_coord,
    upper_left_indices = c(4, 3)
  )

  expect_equal(observed, original_coord)
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
  expect_error(
    crop_points(data.frame(x = 0.5, y_down = 0.5)),
    "Supply `upper_left_indices`"
  )
  expect_error(
    uncrop_points(
      data.frame(x = 0.5, y_down = 0.5),
      upper_left_indices = c(1, 1),
      snapped_zenith_coord = c(2, 2),
      horizon_radius = 1
    ),
    "not both",
    fixed = TRUE
  )
})
