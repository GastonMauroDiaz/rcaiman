histogram_fixture <- function(counts, lower = 0, width = 1) {
  breaks <- lower + (0:length(counts)) * width
  list(breaks = breaks, mids = lower + (seq_along(counts) - .5) * width,
       counts = as.numeric(counts))
}

test_that("histograms accumulate weights without changing observation values", {
  h <- build_histogram(c(0, 1, 2, 3), c(1, 1, 2, 4), n_bins = 3)
  expect_s3_class(h, "histogram")
  expect_equal(h$breaks, 0:3)
  expect_equal(h$mids, c(.5, 1.5, 2.5))
  expect_equal(h$counts, c(1, 1, 6))
  expect_equal(sum(h$density * diff(h$breaks)), 1)
  repeated <- build_histogram(rep(0:3, c(1, 1, 2, 4)), n_bins = 3)
  expect_equal(h$counts, repeated$counts)
  scaled <- build_histogram(0:3, c(1, 1, 2, 4) * .13, 3)
  expect_equal(scaled$counts, h$counts * .13)
  expect_equal(scaled$density, h$density)
  missing <- build_histogram(c(0:3, NA, -100, 100), c(1, 1, 2, 4, 1, 0, NA), 3)
  expect_equal(missing$counts, h$counts)
  expect_equal(missing$breaks, h$breaks)
  for (x in list(matrix(0:3, ncol = 1), data.frame(x = 0:3))) {
    expect_equal(build_histogram(x, data.frame(w = c(1, 1, 2, 4)), 3)$counts, h$counts)
  }
})

test_that("histograms reject undefined ranges, weights, and bin counts", {
  expect_error(build_histogram(c(1, 1)), "distinct")
  expect_error(build_histogram(c(NA, 1)), "distinct")
  expect_error(build_histogram(0:3, rep(0, 4)), "distinct")
  expect_error(build_histogram(0:3, c(1, -1, 1, 1)), "non-negative")
  expect_error(build_histogram(0:3, 1), "length")
  expect_error(build_histogram(c(0, Inf)), "Inf")
  expect_error(build_histogram(matrix(1:4, ncol = 2)), "single column")
  expect_error(build_histogram(c(-1e308, 1e308)), "range")
  for (n in list(1, 2.5, NA_real_, Inf)) expect_error(build_histogram(0:3, n_bins = n))
  expect_error(build_histogram(c(1, 1 + .Machine$double.eps), n_bins = 256), "boundaries")
})

test_that("smoothing retains the grid and mass and handles wide kernels", {
  h <- histogram_fixture(c(10, 1, 0, 0, 0, 1))
  original <- h
  for (sigma in c(.5, 2, 100)) {
    hs <- smooth_histogram(h, sigma)
    expect_identical(h, original)
    expect_identical(hs$breaks, h$breaks)
    expect_identical(hs$mids, h$mids)
    expect_equal(sum(hs$counts), sum(h$counts))
    expect_true(all(is.finite(hs$counts) & hs$counts >= 0))
    expect_equal(sum(hs$density * diff(hs$breaks)), 1)
    reflected <- h
    reflected$counts <- rev(h$counts)
    expect_equal(smooth_histogram(reflected, sigma)$counts, rev(hs$counts))
  }
  expect_identical(smooth_histogram(h, NULL), h)
  expect_length(smooth_histogram(histogram_fixture(c(2, 1)))$counts, 2)
  expect_error(smooth_histogram(h, 0))
})

test_that("local peaks have explicit plateau and endpoint conventions", {
  expect_identical(find_histogram_peaks(histogram_fixture(c(0, 5, 5, 0, 4, 0))), c(2L, 5L))
  expect_identical(find_histogram_peaks(histogram_fixture(c(3, 1, 2, 2, 2, 0, 4))), c(1L, 4L, 7L))
  expect_identical(find_histogram_peaks(histogram_fixture(rep(1, 10))), integer(0))
  expect_identical(find_histogram_peaks(histogram_fixture(c(0, 0, 3, 0, 0))), 3L)
  expect_identical(find_histogram_peaks(histogram_fixture(c(3, 2, 1))), 1L)
})

test_that("Macfarlane windows report bin indices and require 256 bins", {
  counts <- rep(0, 256)
  counts[c(31, 221)] <- c(100, 80)
  h <- histogram_fixture(counts)
  expect_identical(find_macfarlane_peaks(h), c(31L, 221L))
  # Peaks close to the inner limits require another 25-position search step.
  counts[] <- 0
  counts[c(53, 204)] <- 10
  expect_identical(find_macfarlane_peaks(histogram_fixture(counts)), c(53L, 204L))
  expect_error(find_macfarlane_peaks(histogram_fixture(1:32)), "256 bins")
  counts[] <- 0
  counts[128] <- 1
  expect_error(find_macfarlane_peaks(histogram_fixture(counts)), "20 bins")
  counts[] <- 0
  counts[256] <- 1
  expect_error(find_macfarlane_peaks(histogram_fixture(counts)), "suitable peak")
})

test_that("Rosin finds the geometric corner and respects reflection and units", {
  for (n in c(32, 128, 256)) {
    positions <- 0:(n - 1)
    counts <- 1000 * exp(-positions / (n / 8))
    h <- histogram_fixture(counts)
    # Perpendicular distance to one line has the same maximizer as
    # vertical distance. This oracle uses the frequencies without rescaling.
    expected <- h$mids[which.max(counts[1] * (n - positions) / n - counts)]
    expect_equal(thr_rosin(h, 1), expected)
    reflected <- histogram_fixture(rev(counts), lower = -n)
    expect_equal(thr_rosin(reflected, n, "left"), -expected)
    affine <- histogram_fixture(counts, lower = -1, width = 2 / n)
    expect_equal(thr_rosin(affine, 1), -1 + expected * 2 / n)
    h$counts <- counts * .037
    expect_equal(thr_rosin(h, 1), expected)
  }
  # An explicitly supplied peak, not necessarily the global maximum.
  h <- histogram_fixture(c(50, 1, 10, 4, 1, .5, .1, 0))
  expect_true(thr_rosin(h, 3) > h$mids[3])
  expect_error(thr_rosin(h, 8), "positive frequency")
  expect_error(thr_rosin(h, 1, "left"), "No corner")
  for (peak in list(0, 9, 2.5, NA_real_, c(1, 3))) expect_error(thr_rosin(h, peak))
  expect_error(thr_rosin(h, 1, "auto"))
})

test_that("two-corner computes all thresholds from the supplied histogram", {
  counts <- rep(1, 20)
  counts[c(3, 18)] <- 20
  h <- histogram_fixture(counts)
  curve <- h
  h$counts[11] <- 0
  result <- thr_two_corner(h, c(3, 18))
  expect_equal(result$lc, h$mids[4])
  expect_equal(result$uc, h$mids[17])
  expect_equal(result$ts, h$mids[11])
  expect_equal(result$tm, mean(c(result$lc, result$uc)))
  expect_equal(result$tl, result$lc + .25 * (result$uc - result$lc))
  expect_equal(result$th, result$lc + .75 * (result$uc - result$lc))
  # The supplied curve determines ts as well as the corners.
  # This flat valley selects its first bin strictly between the corners.
  expect_equal(thr_two_corner(curve, c(3, 18))$ts, curve$mids[5])
  h$counts[c(3, 18)] <- 0
  expect_error(thr_two_corner(h, c(3, 18)), "positive frequency")
  expect_error(thr_two_corner(curve, c(18, 3)), "increasing order")
  expect_error(thr_two_corner(curve, c(3, 3)), "distinct")
  expect_warning(result <- thr_two_corner(histogram_fixture(c(10, 1, 1, 10)), c(1, 4)), "strictly between")
  expect_true(is.na(result$ts))
  expect_true(is.finite(result$tm))
})

test_that("histogram consumers validate a plain list without requiring a class", {
  h <- histogram_fixture(c(5, 2, 1, 0))
  expect_true(is.finite(thr_rosin(h, 1)))
  for (f in list(smooth_histogram, find_histogram_peaks, function(h) thr_rosin(h, 1))) {
    expect_error(f(list(counts = 1:3)), "breaks")
    invalid <- h
    invalid$breaks[2] <- .5
    expect_error(f(invalid), "equally spaced")
    invalid <- h
    invalid$counts[] <- 0
    expect_error(f(invalid), "positive total")
    invalid$counts[1] <- -1
    expect_error(f(invalid), "non-negative")
  }
})

test_that("diagnostics preserve thresholds and graphical margins", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) })
  graphics::par(mfrow = c(1, 2))
  before <- graphics::par(c("mar", "mfrow"))
  h <- histogram_fixture(exp(-(0:255) / 30))
  expect_equal(thr_rosin(h, 1, diagnose = TRUE), thr_rosin(h, 1))
  h$counts <- rev(h$counts)
  expect_equal(thr_rosin(h, 256, "left", diagnose = TRUE), thr_rosin(h, 256, "left"))
  x <- c(qnorm(ppoints(2000), 30, 7), qnorm(ppoints(1800), 210, 9),
         seq(0, 255, length.out = 200))
  h <- build_histogram(x)
  peaks <- find_macfarlane_peaks(h)
  hs <- smooth_histogram(h)
  for (slope in c(FALSE, TRUE)) {
    expect_equal(thr_two_corner(hs, peaks, slope, diagnose = TRUE),
                 thr_two_corner(hs, peaks, slope))
  }
  expect_identical(graphics::par(c("mar", "mfrow")), before)
})
