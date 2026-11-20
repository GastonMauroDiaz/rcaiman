test_that("outlier thresholds use the default R MAD scaling in every tail", {
  reference <- c(0, 1, 2, 3, 4, 8)
  # Raw MAD would flag 8 at laxity 3; scaled MAD retains it.
  expect_false(any(is_outlier(reference, laxity = 3)))
  threshold <- 3 * stats::mad(reference)
  center <- stats::median(reference)
  x <- center + c(-threshold - .01, -threshold, 0, threshold, threshold + .01)
  expect_identical(is_outlier(x, reference, laxity = 3),
                   c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(is_outlier(x, reference, laxity = 3, cutoff_side = "left"),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(is_outlier(x, reference, laxity = 3, cutoff_side = "right"),
                   c(FALSE, FALSE, FALSE, FALSE, TRUE))
})
