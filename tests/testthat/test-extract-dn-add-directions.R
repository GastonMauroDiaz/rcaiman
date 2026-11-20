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


test_that("extraction functions reject continuous image coordinates", {
  r <- terra::rast(nrows = 3, ncols = 3)
  names(r) <- "dn"
  sampling_coord <- data.frame(x = 0.5, y_down = 0.5)

  expect_error(
    extract_dn(r, sampling_coord, use_window = FALSE),
    "required columns: col, row",
    fixed = TRUE
  )


})
