regional_threshold_fixture <- function() {
  r <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 5))
  r[] <- c(0, 4, 5, 9, 9, 5, 4, 0, 100, 100)
  seg <- terra::rast(r)
  seg[] <- c(rep(2, 4), rep(9, 4), 0, NA)
  w <- terra::rast(r)
  w[] <- c(10, 1, 1, 1, 1, 1, 1, 10, NA, NA)
  list(r = r, seg = seg, w = w)
}

test_that("regional IsoData passes matching weights and retains the binary contract", {
  f <- regional_threshold_fixture()
  result <- binarize_by_region(f$r, f$seg, "thr_isodata", f$w)
  expected <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  expect_identical(as.vector(result[]), expected)
  unweighted <- binarize_by_region(f$r, f$seg, "thr_isodata")
  expect_identical(as.vector(unweighted[]),
    c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(terra::values(binarize_by_region(f$r, f$seg, "thr_isodata", f$w * 17)),
               terra::values(result))
  f$w[] <- 1
  expect_equal(terra::values(binarize_by_region(f$r, f$seg, "thr_isodata", f$w)),
               terra::values(unweighted))
  expect_true(terra::compareGeom(result, f$r))
  expect_false(anyNA(result[]))
})

test_that("zero and missing weights exclude estimation without masking classification", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 4))
  r[] <- c(0, 100, 4, 9)
  seg <- terra::rast(r)
  seg[] <- TRUE
  w <- terra::rast(r)
  w[] <- c(1, NA, 0, 1)
  result <- binarize_by_region(r, seg, "thr_isodata", w)
  expect_identical(as.vector(result[]), c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(terra::values(binarize_by_region(r, as.numeric(seg), "thr_isodata", w)),
               terra::values(result))
})

test_that("insufficient effective observations identify the failing segment", {
  f <- regional_threshold_fixture()
  f$w[5:8] <- c(1, 0, NA, 0)
  expect_error(binarize_by_region(f$r, f$seg, "thr_isodata", f$w),
               "segment 9: At least two observations")
  f$w[1:4] <- 0
  expect_error(binarize_by_region(f$r, f$seg, "thr_isodata", f$w), "segment 2")
  f <- regional_threshold_fixture()
  f$seg[1] <- 7
  expect_error(binarize_by_region(f$r, f$seg, "thr_isodata"), "segment 7")
})

test_that("constant and empty regions still produce logical values without NA", {
  f <- regional_threshold_fixture()
  f$seg[f$seg == 9] <- 2
  f$r[] <- 4
  expect_warning(result <- binarize_by_region(f$r, f$seg, "thr_isodata", f$w),
                  "All retained values are equal")
  expect_identical(as.vector(result[]), rep(FALSE, 10))
  f$seg[] <- 0
  result <- binarize_by_region(f$r, f$seg, "thr_isodata", f$w)
  expect_identical(as.vector(result[]), rep(FALSE, 10))
})

test_that("unsupported methods and invalid weight rasters fail explicitly", {
  f <- regional_threshold_fixture()
  for (method in c("thr_nbhu", "IsoData", "Huang2")) {
    expect_error(binarize_by_region(f$r, f$seg, method, f$w), "supported only")
  }
  for (value in c(-1, Inf)) {
    f$w[1] <- value
    expect_error(binarize_by_region(f$r, f$seg, "thr_isodata", f$w))
  }
  f <- regional_threshold_fixture()
  expect_error(binarize_by_region(f$r, f$seg, "thr_isodata", c(f$w, f$w)), "single layer")
  terra::ext(f$w) <- terra::ext(0, 50, 0, 20)
  expect_error(binarize_by_region(f$r, f$seg, "thr_isodata", f$w), "same geometry")
})

test_that("two-corner failures use the mean of all successful segment thresholds", {
  x <- c(qnorm(ppoints(2000), 30, 7), qnorm(ppoints(1800), 210, 9),
         seq(0, 255, length.out = 200))
  y <- rep(x, 2) + 200
  estimate <- function(x) {
    h <- build_histogram(x)
    thr_two_corner(h, find_macfarlane_peaks(h), slope_reduction = TRUE)$tm
  }
  thresholds <- c(estimate(x), estimate(y))
  fallback <- mean(thresholds)
  values <- c(150, 150, x, rep(250, 3), y, 999, 999)
  labels <- c(7, 7, rep(2, length(x)), rep(9, 3), rep(12, length(y)), 0, NA)
  r <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = length(values)))
  r[] <- values
  seg <- terra::rast(r)
  seg[] <- labels
  w <- terra::rast(r)
  w[] <- ifelse(labels == 12, 10, 1)
  expected <- c(rep(150 > fallback, 2), x > thresholds[1],
                rep(250 > fallback, 3), y > thresholds[2], FALSE, FALSE)
  # The constant segments fail, including the first segment encountered.
  # Unequal sizes and total weights must not weight the fallback mean.
  expect_true(150 < fallback && 250 > fallback)
  testthat::local_mocked_bindings(thr_isodata = function(...) stop("Unexpected IsoData call"))
  expect_message(result <- binarize_by_region(r, seg, "thr_two_corner", w), "mean threshold.*2 segment")
  expect_identical(as.vector(result[]), expected)
  expect_false(anyNA(result[]))
  expect_true(terra::compareGeom(result, r))
  r[] <- rev(values)
  seg[] <- rev(labels)
  w[] <- rev(as.vector(w[]))
  expect_message(reversed <- binarize_by_region(r, seg, "thr_two_corner", w), "mean threshold")
  expect_identical(as.vector(reversed[]), rev(expected))
})

test_that("two-corner has no mean fallback when every segment fails", {
  f <- regional_threshold_fixture()
  f$r[] <- 4
  expect_error(binarize_by_region(f$r, f$seg, "thr_two_corner"), "every segment")
  f$w[] <- 0
  expect_error(binarize_by_region(f$r, f$seg, "thr_two_corner", f$w), "no successful thresholds")
  f$seg[] <- 0
  expect_identical(as.vector(binarize_by_region(f$r, f$seg, "thr_two_corner")[]), rep(FALSE, 10))
})

test_that("regional two-corner uses the weighted histogram pipeline", {
  x <- c(qnorm(ppoints(2000), 30, 7), qnorm(ppoints(1800), 210, 9),
         seq(0, 255, length.out = 200))
  weights <- seq(.1, 3, length.out = length(x))
  h <- build_histogram(x, weights, n_bins = 256)
  threshold <- thr_two_corner(h, find_macfarlane_peaks(h), slope_reduction = TRUE)$tm
  r <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = length(x)))
  r[] <- x
  seg <- terra::rast(r)
  seg[] <- 1
  w <- terra::rast(r)
  w[] <- weights
  testthat::local_mocked_bindings(thr_isodata = function(...) stop("Unexpected fallback"))
  result <- binarize_by_region(r, seg, "thr_two_corner", w)
  expect_identical(as.vector(result[]), x > threshold)
})
