gap_test_raster <- function(x) {
  r <- set_rcaiman_geometry(terra::rast(nrows = nrow(x), ncols = ncol(x)))
  terra::values(r) <- as.vector(t(x))
  names(r) <- "sky"
  r
}

test_that("opening removes isolated gaps and clips rectangle corners", {
  x <- matrix(FALSE, 9, 13)
  x[2:7, 7:11] <- TRUE
  x[3, 2] <- TRUE
  expected <- x
  expected[3, 2] <- FALSE
  expected[cbind(c(2, 2, 7, 7), c(7, 11, 7, 11))] <- FALSE
  r <- gap_test_raster(x)
  result <- rem_small_gaps(r, 3)
  expect_identical(unname(terra::as.matrix(result, wide = TRUE)), expected)
  expect_true(terra::compareGeom(result, r))
  expect_identical(names(result), names(r))
  expect_false(anyNA(terra::values(result)))
  expect_false(any(terra::values(result & !r)))
  expect_equal(terra::values(rem_small_gaps(result, 3)), terra::values(result))
})

test_that("diamond-shaped gaps fitting the structuring element survive", {
  x <- matrix(FALSE, 9, 11)
  for (i in 1:9) for (j in 1:11) x[i, j] <- abs(i - 5) + abs(j - 6) <= 2
  r <- gap_test_raster(x)
  expect_identical(unname(terra::as.matrix(rem_small_gaps(r, 5), wide = TRUE)), x)
  expect_false(any(terra::values(rem_small_gaps(r, 7))))
})

test_that("neutral boundaries preserve constants even on tiny grids", {
  # Both constants, degenerate row/column grids, and an ordinary rectangle.
  cases <- list(c(1, 1, FALSE), c(1, 1, TRUE), c(1, 5, TRUE),
                c(5, 1, FALSE), c(7, 11, TRUE))
  for (case in cases) {
    r <- gap_test_raster(matrix(as.logical(case[3]), case[1], case[2]))
    expect_identical(terra::values(rem_small_gaps(r, 9)), terra::values(r))
  }
  x <- matrix(TRUE, 3, 4); x[1, 1] <- FALSE
  result <- rem_small_gaps(gap_test_raster(x), 3)
  expect_identical(unname(terra::as.matrix(result, wide = TRUE)), x)
})

test_that("size and binary contracts are explicit", {
  r <- gap_test_raster(matrix(TRUE, 3, 5))
  expect_identical(rem_small_gaps(r, 1), r)
  for (size in list(0, -1, 2, 2.5, NA_real_, Inf, c(3, 5))) {
    expect_error(rem_small_gaps(r, size), "size")
  }
  r[2, 2] <- NA
  expect_error(rem_small_gaps(r), "NA-free")
  expect_error(rem_small_gaps(r * 1), "logical")
})
