direction_windows <- function() {
  r <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 5))
  terra::values(r) <- 1:5
  names(r) <- "dn"
  g <- c(r, r)
  names(g) <- c("off_axis", "internal_azimuth")
  terra::values(g) <- cbind(seq(0, 40, 10), 0)
  queries <- data.frame(off_axis = c(10, 30) * pi / 180, internal_azimuth = 0)
  list(r = r, g = g, queries = queries, m = !is.na(r))
}

test_that("overlapping windows combine estimates on the original grid", {
  x <- direction_windows()
  before <- terra::values(x$g)
  out <- apply_by_direction(x$r, x$g, x$queries, x$m, fov = 42, fun = function(v, ...) mean(v))
  expect_identical(names(out), c("v1", "n"))
  expect_equal(terra::values(out$v1, mat = FALSE), c(2.5, 3, 3, 3, 3.5))
  expect_equal(terra::values(out$n, mat = FALSE), c(1, 2, 2, 2, 1))
  expect_true(terra::compareGeom(out, x$r))
  expect_identical(terra::values(x$g), before)
  expect_identical(names(x$m), "dn")
})

test_that("queries are continuous and need not fall inside coverage", {
  x <- direction_windows()
  query <- data.frame(off_axis = 15 * pi / 180, internal_azimuth = 0)
  out <- apply_by_direction(x$r, x$g, query, x$m, fov = 12, fun = function(v, ...) mean(v))
  expect_equal(terra::values(out$v1, mat = FALSE), c(NA, 2.5, 2.5, NA, NA))
  query$off_axis <- 51 * pi / 180
  out <- apply_by_direction(x$r, x$g, query, x$m, fov = 24, fun = function(v, ...) mean(v))
  expect_equal(terra::values(out$v1, mat = FALSE), c(NA, NA, NA, NA, 5))
  expect_equal(terra::values(out$n, mat = FALSE), c(0, 0, 0, 0, 1))
})

test_that("masks and missing data exclude samples and output cells", {
  x <- direction_windows()
  x$m[2] <- FALSE
  x$r[3] <- NA_real_
  x$g[4] <- NA_real_
  query <- data.frame(off_axis = 20 * pi / 180, internal_azimuth = 0)
  out <- apply_by_direction(x$r, x$g, query, x$m, fov = 90, fun = function(v, ...) mean(v))
  expect_equal(terra::values(out$v1, mat = FALSE), c(3, NA, NA, NA, 3))
  expect_equal(terra::values(out$n, mat = FALSE), c(1, 0, 0, 0, 1))
})

test_that("local images preserve spatial arrangement and their valid mask", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 3, ncols = 4))
  terra::values(r) <- 1:12
  names(r) <- "dn"
  g <- c(r * 0 + 20, r * 0)
  names(g) <- c("off_axis", "internal_azimuth")
  m <- r > 4 & r != 6 & r != 8
  r[10] <- NA_real_
  queries <- data.frame(off_axis = 0, internal_azimuth = NA_real_)
  called <- FALSE
  local_fun <- function(local, mask, weights) {
    called <<- TRUE
    expect_equal(dim(local), c(2, 4, 1))
    expect_identical(names(local), "dn")
    expect_equal(unname(as.vector(terra::ext(local))), c(0, 4, 0, 2))
    expect_identical(terra::values(mask, mat = FALSE),
                     c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE))
    expect_equal(terra::values(local, mat = FALSE), c(5, NA, 7, NA, 9, NA, 11, 12))
    expect_true(terra::compareGeom(local, weights))
    expect_equal(terra::values(weights, mat = FALSE), c(1, NA, 1, NA, 1, NA, 1, 1))
    mean(terra::values(local)[terra::values(mask)])
  }
  out <- apply_by_direction(r, g, queries, m, fov = 80,
                            method = "single_layer_image", fun = local_fun)
  expect_true(called)
  vector <- apply_by_direction(r, g, queries, m, fov = 80, fun = function(v, ...) mean(v))
  expect_equal(terra::values(out), terra::values(vector))
})

test_that("values and unnormalized weights reach the same window in cell order", {
  x <- direction_windows()
  w <- terra::setValues(terra::rast(x$r), c(1, 3, 9, 2, 5))
  seen <- list()
  fun <- function(v, weights) {
    seen[[length(seen) + 1L]] <<- list(v = v, w = weights)
    stats::weighted.mean(v, weights)
  }
  out <- apply_by_direction(x$r, x$g, x$queries, x$m, fun, w = w, fov = 42)
  expect_equal(seen, list(list(v = 1:4, w = c(1, 3, 9, 2)),
                          list(v = 2:5, w = c(3, 9, 2, 5))))
  first <- weighted.mean(1:4, c(1, 3, 9, 2))
  second <- weighted.mean(2:5, c(3, 9, 2, 5))
  expect_equal(terra::values(out$v1, mat = FALSE),
               c(first, rep(mean(c(first, second)), 3), second))
  expect_equal(terra::values(out$n, mat = FALSE), c(1, 2, 2, 2, 1))
  scaled <- apply_by_direction(x$r, x$g, x$queries, x$m, weighted.mean,
                               w = w * 100, fov = 42)
  expect_equal(terra::values(scaled), terra::values(out))
  unweighted <- apply_by_direction(x$r, x$g, x$queries, x$m,
                                    fun = function(v, w) { expect_equal(w, rep(1, length(v))); mean(v) },
                                    fov = 42)
  expect_equal(terra::values(unweighted$v1, mat = FALSE), c(2.5, 3, 3, 3, 3.5))
})

test_that("weight rasters preserve local alignment and missing weights exclude pairs", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 3, ncols = 4))
  terra::values(r) <- 1:12
  g <- c(r * 0 + 20, r * 0)
  names(g) <- c("off_axis", "internal_azimuth")
  w <- r^2
  w[6] <- NA_real_
  w[7] <- 0
  r[10] <- NA_real_
  m <- r > 4
  m[10] <- TRUE
  q <- data.frame(off_axis = 0, internal_azimuth = NA_real_)
  fun <- function(r, m, w) {
    expect_equal(dim(r), c(2, 4, 1))
    expect_true(terra::compareGeom(r, w))
    expect_identical(terra::values(m, mat = FALSE), !is.na(terra::values(w, mat = FALSE)))
    expect_equal(terra::values(r, mat = FALSE), c(5, NA, 7, 8, 9, NA, 11, 12))
    expect_equal(terra::values(w, mat = FALSE), c(25, NA, 0, 64, 81, NA, 121, 144))
    selected <- terra::values(m, mat = FALSE)
    weighted.mean(terra::values(r, mat = FALSE)[selected],
                  terra::values(w, mat = FALSE)[selected])
  }
  image <- apply_by_direction(r, g, q, m, fun, w = w,
                              method = "single_layer_image", fov = 80)
  vector <- apply_by_direction(r, g, q, m, weighted.mean, w = w, fov = 80)
  expect_equal(terra::values(image), terra::values(vector))
  expect_equal(terra::values(image$n, mat = FALSE), c(rep(0, 4), 1, 0, 1, 1, 1, 0, 1, 1))
})

test_that("supplying weights does not change the combination of window estimates", {
  x <- direction_windows()
  fun <- function(v, ...) mean(v)
  baseline <- apply_by_direction(x$r, x$g, x$queries, x$m, fun, fov = 42)
  supplied <- apply_by_direction(x$r, x$g, x$queries, x$m, fun,
                                 w = x$r^3, fov = 42)
  expect_equal(terra::values(supplied), terra::values(baseline))
  # The two-corner-style wrapper controls errors while accepting weights.
  unavailable <- function(v, ...) tryCatch(stop("cannot fit"), error = function(e) NA)
  out <- apply_by_direction(x$r, x$g, x$queries, x$m, unavailable, w = x$r)
  expect_true(all(is.na(terra::values(out$v1))))
  expect_true(all(terra::values(out$n) == 0))
})

test_that("weights follow the raster contract and the estimator is explicit", {
  x <- direction_windows()
  expect_error(apply_by_direction(x$r, x$g, x$queries, x$m), "fun")
  w <- x$r
  w[1] <- -1
  expect_error(apply_by_direction(x$r, x$g, x$queries, x$m, weighted.mean, w), "non-negative")
  w <- terra::setValues(terra::rast(x$r), c(Inf, 2, 3, 4, 5))
  expect_error(apply_by_direction(x$r, x$g, x$queries, x$m, weighted.mean, w), "Inf")
  expect_error(apply_by_direction(x$r, x$g, x$queries, x$m, weighted.mean, x$g), "single layer")
  w <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 3))
  w[] <- 1
  expect_error(apply_by_direction(x$r, x$g, x$queries, x$m, weighted.mean, w), "same geometry")
  w <- x$r * 0
  called <- FALSE
  out <- apply_by_direction(x$r, x$g, x$queries, x$m,
                            function(v, w) { called <<- TRUE; expect_true(all(w == 0)); NA },
                            w = w)
  expect_true(called)
  expect_true(all(terra::values(out$n) == 0))
})

test_that("each estimator has its own contributions and count", {
  x <- direction_windows()
  fun <- function(v, ...) c(mean = mean(v), other = if (min(v) == 1) NA_real_ else 9)
  out <- apply_by_direction(x$r, x$g, x$queries, x$m, fov = 42, fun = fun)
  expect_identical(names(out), c("mean", "other", "n_mean", "n_other"))
  expect_equal(terra::values(out$other, mat = FALSE), c(NA, 9, 9, 9, 9))
  expect_equal(terra::values(out$n_mean, mat = FALSE), c(1, 2, 2, 2, 1))
  expect_equal(terra::values(out$n_other, mat = FALSE), c(0, 1, 1, 1, 1))
  # Named results can change order without changing their meaning.
  fun <- function(v, ...) if (min(v) == 1) c(a = 1, b = 2) else c(b = 2, a = 1)
  out <- apply_by_direction(x$r, x$g, x$queries, x$m, fov = 42, fun = fun)
  expect_equal(terra::values(out$a, mat = FALSE), rep(1, 5))
  expect_equal(terra::values(out$b, mat = FALSE), rep(2, 5))
})

test_that("MAD filtering is retained and counts precede rejection", {
  x <- direction_windows()
  queries <- x$queries[rep(1, 5), ]
  index <- 0L
  fun <- function(v, ...) { index <<- index + 1L; c(1, 2, 3, 4, 100)[index] }
  out <- apply_by_direction(x$r, x$g, queries, x$m, fov = 360, fun = fun)
  expect_equal(terra::values(out$v1, mat = FALSE), rep(2.5, 5))
  expect_equal(terra::values(out$n, mat = FALSE), rep(5, 5))
  index <- 0L
  out <- apply_by_direction(x$r, x$g, queries[1:2, ], x$m, fov = 360,
                            laxity = 0.1, fun = fun)
  expect_true(all(is.na(terra::values(out$v1))))
  expect_equal(terra::values(out$n, mat = FALSE), rep(2, 5))
})

test_that("empty windows and unavailable estimates have simple outputs", {
  x <- direction_windows()
  fail <- function(v, ...) stop("should not be called")
  for (queries in list(x$queries[FALSE, ],
                       data.frame(off_axis = pi, internal_azimuth = 0),
                       data.frame(off_axis = NA_real_, internal_azimuth = 0))) {
    out <- apply_by_direction(x$r, x$g, queries, x$m, fov = 20, fun = fail)
    expect_identical(names(out), c("v1", "n"))
    expect_true(all(is.na(terra::values(out$v1))))
    expect_true(all(terra::values(out$n) == 0))
  }
  out <- apply_by_direction(x$r, x$g, x$queries, x$m & FALSE, fun = fail)
  expect_true(all(terra::values(out$n) == 0))
  for (fun in list(function(v, ...) NA, function(v, ...) NULL, function(v, ...) numeric(0))) {
    out <- apply_by_direction(x$r, x$g, x$queries, x$m, fun = fun)
    expect_true(all(is.na(terra::values(out$v1))))
    expect_true(all(terra::values(out$n) == 0))
  }
  out <- apply_by_direction(x$r, x$g, x$queries, x$m, fun = function(v, ...) c(a = NA, b = NA))
  expect_identical(names(out), c("a", "b", "n_a", "n_b"))
  expect_true(all(is.na(terra::values(out[[1:2]]))))
  expect_true(all(terra::values(out[[3:4]]) == 0))
})

test_that("invalid results report which query failed", {
  x <- direction_windows()
  call <- function(fun, ...) apply_by_direction(x$r, x$g, x$queries, x$m,
                                               fov = 42, fun = fun, ...)
  expect_error(call(function(v, ...) stop("cannot fit")), "row 1: cannot fit")
  expect_error(call(function(r, m, w) stop("image failed"), method = "single_layer_image"),
                "row 1: image failed")
  expect_error(call(function(v, ...) "oops"), "numeric vector.*row 1")
  expect_error(call(function(v, ...) Inf), "finite or NA.*row 1")
  expect_error(call(function(v, ...) c(a = 1, a = 2)), "unique.*row 1")
  expect_error(call(function(v, ...) c(n = 2)), "collide")
  expect_error(call(function(v, ...) c(a = 1, n_a = 2)), "collide")
  expect_error(call(function(v, ...) if (min(v) == 1) c(a = 1) else c(b = 2)),
                "inconsistent.*row 2")
  expect_error(call(2), "must be a function")
  expect_error(apply_by_direction(x$r, x$g, x$queries, x$m, fov = 361), "360")
  bad <- x$queries
  bad$off_axis <- bad$off_axis * 180 / pi
  expect_error(apply_by_direction(x$r, x$g, bad, x$m, fun = function(v, ...) mean(v)), "radians")
  bad <- x$g
  names(bad) <- c("zenith", "azimuth")
  expect_error(apply_by_direction(x$r, bad, x$queries, x$m), "`g`")
})
