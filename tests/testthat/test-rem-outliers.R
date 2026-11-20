outlier_trend_samples <- function() {
  x <- sky_grid_centers(30)
  x$id <- seq_len(nrow(x))
  x$Blue <- with(x, 100 + 200 * zenith^2 + 5 * sin(azimuth) -
                  3 * cos(azimuth) + .2 * sin(id))
  x$Blue[1] <- x$Blue[1] + 50
  x
}

test_that("query identity is excluded even with coincident directions", {
  x <- data.frame(zenith = rep(.5, 5), azimuth = rep(.3, 5),
                  dn = c(0, 1, 2, 100, 200), output = letters[1:5])
  rownames(x) <- paste0("sample_", 1:5)
  got <- rem_outliers(x, "dn", k = 3, search_radius = NULL)
  expect_identical(got, x[1:3, ])
  expect_identical(rem_outliers(x, 3, k = 3, search_radius = NULL), got)
  # Exactly three OTHER observations are required, rather than two plus self.
  expect_identical(rem_outliers(x[1:3, ], "dn", k = 3), x[1:3, ])
  expect_identical(rem_outliers(x[1:4, ], "dn", k = 4), x[1:4, ])
  expect_identical(rem_outliers(x[1:4, ], "dn", k = 3), x[1:3, ])
})

test_that("angular caps retain isolated or undefined queries", {
  x <- data.frame(zenith = c(rep(pi/2, 5), NA),
                  azimuth = c(0:3, 120, 0) * pi/180,
                  dn = c(0, 1, 2, 100, 1000, 2000))
  expect_equal(rem_outliers(x, "dn", k = 3, search_radius = 4), x[-4, ])
  # Azimuth need not be defined for a pole in the geometric search.
  x$zenith <- 0; x$azimuth <- NA_real_
  expect_equal(rem_outliers(x, "dn", k = 3), x[1:3, ])
})

test_that("MAD and cutoff sides agree with explicit leave-one-out comparisons", {
  x <- data.frame(zenith = seq(.1, .9, length.out = 8), azimuth = .5,
                  dn = c(-100, 1, 2, 3, 4, 5, 6, 100))
  for (side in c("both", "left", "right")) {
    keep <- vapply(seq_len(nrow(x)), function(i)
      !is_outlier(x$dn[i], x$dn[-i], cutoff_side = side), logical(1))
    expect_equal(rem_outliers(x, "dn", k = 7, search_radius = NULL, cutoff_side = side),
                 x[keep, ])
  }
  x$dn <- 1
  x$dn[1] <- 100
  # Keep the established no-dispersion rule, even for a high query value.
  expect_equal(rem_outliers(x, "dn", k = 7, search_radius = NULL), x)
})

test_that("detrending uses the selected variable and predicts an excluded query", {
  x <- outlier_trend_samples()
  keep <- vapply(seq_len(nrow(x)), function(i) {
    fit <- lm(Blue ~ zenith + I(zenith^2) + sin(azimuth) + cos(azimuth), data = x[-i, ])
    !is_outlier(x$Blue[i] - as.numeric(predict(fit, newdata = x[i, ])), residuals(fit))
  }, logical(1))
  got <- rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = TRUE)
  expect_equal(got, x[keep, ])
  expect_false(1 %in% got$id)
  raw <- rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL)
  expect_true(1 %in% raw$id)
  expect_equal(rem_outliers(x, 4, k = nrow(x)-1, search_radius = NULL, detrend = TRUE), got)
})

test_that("unsupported local trends stop instead of changing the comparison", {
  x <- outlier_trend_samples()
  expect_error(rem_outliers(x, "Blue", k = 19, search_radius = NULL, detrend = TRUE),
               "at least 20")
  x$zenith <- .5
  expect_error(rem_outliers(x, "Blue", k = 20, search_radius = NULL, detrend = TRUE),
               "row 1: Detrending failed:.*full-rank")
  x <- outlier_trend_samples()
  x$zenith[1] <- 0; x$azimuth[1] <- NA_real_
  expect_error(rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = TRUE),
               "row 1: Detrending failed:.*prediction")
  # The missing external angle can also belong to a neighbor.
  x <- outlier_trend_samples()
  x$zenith[2] <- 0; x$azimuth[2] <- NA_real_
  expect_error(rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = TRUE),
               "row 1: Detrending failed:.*azimuth")
})

test_that("fitting and prediction errors retain their cause and input row", {
  x <- outlier_trend_samples()
  local_mocked_bindings(fit_cone_shaped_model = function(...) stop("test fitting failure"))
  expect_error(rem_outliers(x, "Blue", k = 20, search_radius = NULL, detrend = TRUE),
               "row 1: Detrending failed: test fitting failure")
})

test_that("prediction failure is not replaced by an untransformed comparison", {
  x <- outlier_trend_samples()
  local_mocked_bindings(predict = function(...) stop("test prediction failure"), .package = "stats")
  expect_error(rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = TRUE),
               "row 1: Detrending failed: test prediction failure")
})

test_that("warnings from usable trends remain visible with row context", {
  x <- outlier_trend_samples()
  expected <- rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = TRUE)
  original_fit <- fit_cone_shaped_model
  local_mocked_bindings(fit_cone_shaped_model = function(...) {
    warning("test usable fit warning", call. = FALSE)
    original_fit(...)
  })
  messages <- character()
  got <- withCallingHandlers(
    rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = TRUE),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_identical(got, expected)
  expect_identical(messages, paste0("Detrending at sampling_directions row ",
                                   seq_len(nrow(x)), ": test usable fit warning"))
})

test_that("sparse neighborhoods and empty inputs still retain rows without testing", {
  x <- outlier_trend_samples()[1:5, ]
  local_mocked_bindings(fit_cone_shaped_model = function(...) stop("Should not fit"))
  expect_identical(rem_outliers(x, "Blue", k = 20, detrend = TRUE), x)
  expect_identical(rem_outliers(x[FALSE, ], "Blue", k = 20, detrend = TRUE), x[FALSE, ])
})

test_that("internal directions can select neighbors but trend stays external", {
  x <- outlier_trend_samples()
  x <- internal_directions_from_external(x, sensor_rotation = 17,
                                         azimuth_rotation = 24, elevation = 32.5)
  for (detrend in c(FALSE, TRUE)) {
    external <- rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL, detrend = detrend)
    internal <- rem_outliers(x, "Blue", k = nrow(x)-1, search_radius = NULL,
      detrend = detrend, angle_columns = c("off_axis", "internal_azimuth"))
    expect_equal(internal, external)
  }
  only_internal <- x[c("off_axis", "internal_azimuth", "Blue")]
  expect_no_error(rem_outliers(only_internal, "Blue", k = 3,
    angle_columns = c("off_axis", "internal_azimuth")))
  expect_error(rem_outliers(only_internal, "Blue", detrend = TRUE,
    angle_columns = c("off_axis", "internal_azimuth")), "zenith")
})

test_that("contracts are validated even for empty or sparse tables", {
  x <- outlier_trend_samples()
  expect_identical(rem_outliers(x[FALSE, ], "Blue"), x[FALSE, ])
  expect_identical(rem_outliers(x[1, ], "Blue"), x[1, ])
  expect_error(rem_outliers(x, "missing"), "required columns")
  expect_error(rem_outliers(x, 2.5), "integer")
  expect_error(rem_outliers(x, "Blue", k = 2), "at least 3")
  expect_error(rem_outliers(x, "Blue", detrend = c(TRUE, FALSE)), "length one|length 1")
  expect_error(rem_outliers(x, "Blue", search_radius = 181), "180")
  expect_error(rem_outliers(x, "Blue", angle_columns = c("zenith", "zenith")), "different")
  expect_error(rem_outliers(x[FALSE, ], "Blue", cutoff_side = "invalid"), "cutoff_side")
  expect_error(rem_outliers(x[FALSE, ], "Blue", laxity = 0), "positive")
  for (bad in c(NA, Inf, -Inf)) {
    x$Blue[1] <- bad
    expect_error(rem_outliers(x, "Blue"), "Blue")
  }
})

test_that("parallel search and evaluation match sequential decisions", {
  skip_on_cran()
  skip_if(is.na(parallel::detectCores()) || parallel::detectCores() < 2)
  on.exit(foreach::registerDoSEQ(), add = TRUE)
  x <- outlier_trend_samples()
  # Detrending also performs the MAD comparison, while exercising model
  # fitting and prediction in each worker. Raw comparisons are tested above.
  k <- nrow(x)-1
  expected <- rem_outliers(x, "Blue", k = k, search_radius = NULL, detrend = TRUE)
  actual <- rem_outliers(x, "Blue", k = k, search_radius = NULL, detrend = TRUE,
                         parallel = TRUE, cores = 2)
  expect_identical(actual, expected)
})


test_that("k = 20 supplies twenty training observations for detrending", {
  x <- data.frame(zenith = seq(.15, 1.25, length.out = 21),
                  azimuth = (seq_len(21) * .7) %% (2*pi), id = seq_len(21))
  x$dn <- with(x, 100 + 200 * zenith^2 + 5*sin(azimuth) + .1*sin(id))
  x$dn[1] <- x$dn[1] + 30
  keep <- vapply(seq_len(nrow(x)), function(i) {
    fit <- lm(dn ~ zenith + I(zenith^2) + sin(azimuth) + cos(azimuth), data = x[-i, ])
    !is_outlier(x$dn[i] - as.numeric(predict(fit, newdata = x[i, ])), residuals(fit))
  }, logical(1))
  got <- rem_outliers(x, "dn", k = 20, search_radius = NULL, detrend = TRUE)
  expect_equal(got, x[keep, ])
  expect_false(1 %in% got$id)
  expect_true(1 %in% rem_outliers(x, "dn", k = 20, search_radius = NULL)$id)
})
