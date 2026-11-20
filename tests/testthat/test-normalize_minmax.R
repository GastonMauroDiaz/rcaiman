test_that("normalization uses the declared range and preserves missing values", {
  r <- terra::rast(nrows = 2, ncols = 3, vals = c(0, 64, 128, 192, 255, NA))
  normalized <- normalize_min_max(r, 0, 255)
  expect_equal(terra::values(normalized), terra::values(r) / 255)
})
