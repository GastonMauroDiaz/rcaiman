test_that("IsoData weights group means while keeping the threshold in DN", {
  x <- c(0, 2, 10, 12)
  w <- c(1, 3, 1, 1)
  expect_equal(thr_isodata(x), 6)
  expect_equal(thr_isodata(x, w), 6.25)
  expect_equal(thr_isodata(x, w * 100), 6.25)
  expect_equal(thr_isodata(x, w * 1e-250), 6.25)
  expect_equal(thr_isodata(x, w * 1e300), 6.25)
  expect_equal(thr_isodata(x, rep(7, length(x))), thr_isodata(x))
  expect_equal(thr_isodata(10 + 3 * x, w), 10 + 3 * 6.25)
  expect_equal(thr_isodata(c(-4, -2, 6)), 1.5)
})

test_that("integer weights reproduce repeated observations", {
  x <- c(-3, 0, 2, 8, 10, 20)
  w <- c(1, 3, 2, 2, 5, 1)
  # Independently evaluate the unweighted IsoData iteration on replicated data.
  repeated <- rep(x, w)
  threshold <- mean(repeated)
  repeat {
    updated <- (mean(repeated[repeated <= threshold]) +
                  mean(repeated[repeated > threshold])) / 2
    if (updated == threshold) break
    threshold <- updated
  }
  expect_equal(thr_isodata(x, w), threshold)
  expect_equal(thr_isodata(repeated), threshold)
})

test_that("missing and zero-weight observations are excluded as pairs", {
  x <- c(0, 2, 10, 12)
  w <- c(1, 3, 1, 1)
  expect_equal(thr_isodata(c(NA, x, 1000, -1000), c(100, w, NA, 0)), 6.25)
  expect_equal(thr_isodata(c(x, NA)), 6)
  expect_error(thr_isodata(x, rep(0, 4)), "At least two")
  expect_error(thr_isodata(x, c(1, 0, 0, 0)), "At least two")
  expect_error(thr_isodata(c(NA_real_, NA_real_)), "At least two")
  expect_error(thr_isodata(numeric()), "At least two")
  expect_warning(result <- thr_isodata(c(3, 3, 100), c(1, 2, 0)), "equal")
  expect_equal(result, 3)
})

test_that("single-column inputs work and invalid inputs fail explicitly", {
  x <- c(0, 2, 10, 12)
  w <- c(1, 3, 1, 1)
  expect_equal(thr_isodata(matrix(x), matrix(w)), 6.25)
  expect_equal(thr_isodata(data.frame(dn = x), data.frame(area = w)), 6.25)
  expect_error(thr_isodata(cbind(x, x)), "single column")
  expect_error(thr_isodata(x, cbind(w, w)), "single column")
  expect_error(thr_isodata(as.character(x)), "numeric")
  expect_error(thr_isodata(data.frame(dn = factor(x))), "numeric")
  expect_error(thr_isodata(x, w[-1]), "length")
  expect_error(thr_isodata(x, -w), "non-negative")
  expect_error(thr_isodata(c(x, Inf)), "Inf")
  expect_error(thr_isodata(x, c(1, 1, 1, Inf)), "Inf")
})

test_that("all observations are used without consuming random numbers", {
  set.seed(135)
  before <- .Random.seed
  x <- c(rep(0, 1000), 10, 12)
  expect_equal(thr_isodata(x), 5.5)
  expect_identical(.Random.seed, before)
  expect_equal(thr_isodata(rev(x)), 5.5)
  expect_warning(threshold <- thr_isodata(rep(4, 600)), "equal")
  expect_equal(threshold, 4)
  expect_identical(.Random.seed, before)
})

test_that("area weights pass directly through apply_by_direction to IsoData", {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
                   center_of_vision = c(10, 10), normalization_radius = 10,
                   normalization_angle = 90)
  g <- set_internal_geometry(c(20, 20), geometry)
  r <- g$off_axis + 20
  w <- area_weight_image(g, geometry)
  m <- !is.na(g$off_axis)
  query <- data.frame(off_axis = 0, internal_azimuth = 0)
  expected <- thr_isodata(r[], w[])
  out <- apply_by_direction(r, g, query, m, w = w, fov = 360,
                            fun = thr_isodata)
  selected <- terra::values(m, mat = FALSE)
  expect_equal(terra::values(out$v1, mat = FALSE)[selected],
               rep(expected, sum(selected)))
  expect_true(all(is.na(terra::values(out$v1, mat = FALSE)[!selected])))
  expect_true(all(terra::values(out$n, mat = FALSE)[selected] == 1))
  unweighted <- apply_by_direction(r, g, query, m, fov = 360,
                                   fun = thr_isodata)
  expect_equal(terra::values(unweighted$v1, mat = FALSE)[selected],
               rep(thr_isodata(r[]), sum(selected)))
  expect_gt(abs(expected - thr_isodata(r[])), 0.1)
})
