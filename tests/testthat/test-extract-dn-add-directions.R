test_that("extract_dn() extracts rectangular raster cells in input order", {
  r <- terra::rast(
    nrows = 3,
    ncols = 5,
    xmin = 0,
    xmax = 5,
    ymin = 0,
    ymax = 3
  )
  terra::values(r) <- 1:15
  names(r) <- "dn"

  sampling_indices <- data.frame(
    col = c(1, 1, 2, 5),
    row = c(1, 1, 2, 3),
    id = letters[1:4]
  )

  observed <- extract_dn(r, sampling_indices, use_window = FALSE)

  expect_equal(
    observed,
    data.frame(
      col = sampling_indices$col,
      row = sampling_indices$row,
      id = sampling_indices$id,
      dn = c(1, 1, 7, 15)
    )
  )
})


test_that("extract_dn() preserves the interior 3 x 3 mean", {
  r <- terra::rast(
    nrows = 3,
    ncols = 5,
    xmin = 0,
    xmax = 5,
    ymin = 0,
    ymax = 3
  )
  terra::values(r) <- 1:15
  names(r) <- "dn"

  observed <- extract_dn(
    r,
    data.frame(col = 3, row = 2),
    use_window = TRUE
  )

  expect_equal(observed$dn, 8)
})


test_that("extract_dn() returns NA for windows crossing raster borders", {
  r <- terra::rast(
    nrows = 3,
    ncols = 5,
    xmin = 0,
    xmax = 5,
    ymin = 0,
    ymax = 3
  )
  terra::values(r) <- 1:15
  r <- c(r, r * 10)
  names(r) <- c("dn1", "dn2")

  observed <- extract_dn(
    r,
    data.frame(col = 1, row = 1),
    use_window = TRUE
  )

  expect_true(is.na(observed$dn1))
  expect_true(is.na(observed$dn2))
})


test_that("extract_dn() returns NA for rasters without values", {
  r <- terra::rast(nrows = 3, ncols = 3)
  names(r) <- "dn"
  sampling_indices <- data.frame(col = 2, row = 2)

  observed_cell <- extract_dn(
    r,
    sampling_indices,
    use_window = FALSE
  )
  observed_window <- extract_dn(
    r,
    sampling_indices,
    use_window = TRUE
  )

  expect_true(is.na(observed_cell$dn))
  expect_true(is.na(observed_window$dn))
})


test_that("extraction functions accept zero-row index data", {
  r <- terra::rast(nrows = 3, ncols = 3)
  terra::values(r) <- 1:9
  names(r) <- "dn"
  sampling_indices <- data.frame(
    col = integer(),
    row = integer(),
    id = character()
  )

  observed_cell <- extract_dn(
    r,
    sampling_indices,
    use_window = FALSE
  )
  observed_window <- extract_dn(
    r,
    sampling_indices,
    use_window = TRUE
  )

  expect_identical(names(observed_cell), c("col", "row", "id", "dn"))
  expect_identical(names(observed_window), c("col", "row", "id", "dn"))
  expect_equal(nrow(observed_cell), 0)
  expect_equal(nrow(observed_window), 0)

  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  observed_directions <- add_directions(sampling_indices, z, a)

  expect_identical(
    names(observed_directions),
    c("col", "row", "id", "zen", "azi")
  )
  expect_equal(nrow(observed_directions), 0)
})


test_that("extract_dn() requires integer-like positive indices", {
  r <- terra::rast(nrows = 3, ncols = 3)
  terra::values(r) <- 1:9

  expect_error(
    extract_dn(
      r,
      data.frame(col = 1.5, row = 1),
      use_window = FALSE
    ),
    "sampling_indices\\$col.*integer-like"
  )
  expect_error(
    extract_dn(
      r,
      data.frame(col = 1, row = 0),
      use_window = FALSE
    ),
    "sampling_indices\\$row.*positive"
  )
})


test_that("add_directions() preserves indices and adds angles", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  sampling_indices <- data.frame(
    col = c(5, 6, 8),
    row = c(5, 5, 6),
    id = letters[1:3]
  )

  observed <- add_directions(sampling_indices, z, a)

  expect_equal(observed$col, sampling_indices$col)
  expect_equal(observed$row, sampling_indices$row)
  expect_identical(observed$id, sampling_indices$id)
  expect_equal(
    observed$zen,
    c(0.3141592654, 0.3141592654, 0.9424777961),
    tolerance = 1e-8
  )
  expect_equal(
    observed$azi,
    c(0.7853981634, 5.4977871438, 4.5149934205),
    tolerance = 1e-8
  )
})


test_that("add_directions() preserves additional column names", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  sampling_indices <- data.frame(col = 5, row = 5)
  sampling_indices[["Zenith image note"]] <- "keep zen"
  sampling_indices[["Azimuth image note"]] <- "keep azi"

  observed <- add_directions(sampling_indices, z, a)

  expect_identical(
    names(observed),
    c(
      "col",
      "row",
      "Zenith image note",
      "Azimuth image note",
      "zen",
      "azi"
    )
  )
})
test_that("add_directions() rejects indices outside the angular domain", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)

  expect_error(
    add_directions(data.frame(col = 1, row = 1), z, a),
    "Some sampling points fall outside"
  )
})


test_that("extraction functions reject continuous image coordinates", {
  r <- terra::rast(nrows = 3, ncols = 3)
  names(r) <- "dn"
  sampling_coord <- data.frame(x = 0.5, y_down = 0.5)

  expect_error(
    extract_dn(r, sampling_coord, use_window = FALSE),
    "required columns: col, row",
    fixed = TRUE
  )

  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  expect_error(
    add_directions(sampling_coord, z, a),
    "required columns: col, row",
    fixed = TRUE
  )
})
