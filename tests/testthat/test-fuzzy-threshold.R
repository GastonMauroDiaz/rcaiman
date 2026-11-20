fuzzy_test_raster <- function(x) {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2)
  terra::values(r) <- x
  r
}

test_that("local thresholds move the midpoint and local scales change the slope", {
  r <- fuzzy_test_raster(c(-3, 1, 7, -3, 1, 7))
  threshold <- fuzzy_test_raster(c(-2, 1, 6, -1, 1, 5))
  width <- fuzzy_test_raster(rep(c(1, 2), each = 3) / log(3))
  out <- fuzzy_threshold(r, threshold, width)
  expect_equal(as.vector(out[]), rep(c(.25, .5, .75), 2))
  expect_true(terra::compareGeom(r, out))
  expect_identical(names(out), "membership")
  expect_equal(fuzzy_threshold(r, 1, 2)[], stats::plogis((r[] - 1) / 2),
               ignore_attr = TRUE)
})

test_that("zero scale implements strict thresholding and missing inputs propagate", {
  r <- fuzzy_test_raster(c(0, 1, 2, 1, NA, 2))
  threshold <- fuzzy_test_raster(c(1, 1, 1, 1, 1, NA))
  width <- fuzzy_test_raster(c(0, 0, 0, 1, 0, 0))
  expect_equal(as.vector(fuzzy_threshold(r, threshold, width)[]),
               c(0, 0, 1, .5, NA, NA))
  expect_equal(as.vector(fuzzy_threshold(r, 1, 0)[]), c(0, 0, 1, 0, NA, 1))
  width[1] <- NA
  expect_true(is.na(fuzzy_threshold(r, 1, width)[][1]))
})

test_that("threshold and scale rasters must match the image grid", {
  r <- fuzzy_test_raster(1:6)
  shifted <- terra::shift(r, dx = 1)
  expect_error(fuzzy_threshold(r, shifted, 1), "same geometry")
  expect_error(fuzzy_threshold(r, 1, shifted), "same geometry")
  expect_error(fuzzy_threshold(r, c(r, r), 1), "single layer")
  expect_error(fuzzy_threshold(r, 1, -r), "non-negative")
  expect_error(fuzzy_threshold(r, c(1, 2), 1), "length")
  expect_error(fuzzy_threshold(r, 1, Inf), "Inf")
})
