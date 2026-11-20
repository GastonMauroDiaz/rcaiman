spherical_samples <- function(angle = c(0, 10, 30), dn = c(2, 8, 20)) {
  data.frame(zenith = angle * pi / 180, azimuth = 0, dn = dn)
}
spherical_queries <- function(angle = c(0, 5, 90)) {
  data.frame(zenith = angle * pi / 180, azimuth = 0,
             background = seq_along(angle) * 100)
}

test_that("pure IDW preserves queries and handles exact and absent neighbors", {
  samples <- spherical_samples()
  query <- spherical_queries()
  query$output <- letters[1:3]
  rownames(query) <- c("first", "second", "last")
  out <- interpolate_spherical(samples, query, "dn", k = 2, search_radius = 20)
  expect_identical(out[names(query)], query)
  expect_equal(out$dn, c(2, 5, NA))
  duplicated <- spherical_samples(c(0, 0, 10), c(2, 6, 100))
  out <- interpolate_spherical(duplicated, query, "dn", k = 3)
  expect_equal(out$dn[1], 4)
  for (p in c(.5, 1, 2, 4)) {
    out <- interpolate_spherical(samples, spherical_queries(4), 3L, k = 2, p = p)
    expected <- sum(c(2, 8) / c(4, 6)^p) / sum(1 / c(4, 6)^p)
    expect_equal(out$dn, expected)
  }
})

test_that("fixed-bound blending retains its distance-dependent weights", {
  samples <- spherical_samples(c(0, .5, 1, 5, 20), c(2, 4, 8, 10, 50))
  query <- spherical_queries(0)
  query$background <- 100
  for (p in c(.5, 1, 2, 4)) {
    out <- interpolate_spherical(samples, query, "dn", "background", .1,
                                  k = NULL, p = p, min_dist = 1, search_radius = 20)
    weight <- c(1, 1, 1, (5^(-p) - 20^(-p)) / (1 - 20^(-p)), 0)
    expected <- (sum(weight * samples$dn) + .1 * 100) / (sum(weight) + .1)
    expect_equal(out$dn, expected, tolerance = 1e-10)
  }
  # A coincident observation has finite weight when blending is active.
  out <- interpolate_spherical(samples[1, ], query, "dn", 3L, filling_weight = 1)
  expect_equal(out$dn, 51)
  out <- interpolate_spherical(samples[5, ], query, "dn", "background",
                                search_radius = 20)
  expect_equal(out$dn, 100, tolerance = 1e-10)
  near <- interpolate_spherical(samples[4, ], query, "dn", "background", k = 1)
  far_sample <- samples[4, ]; far_sample$zenith <- 15 * pi / 180
  far <- interpolate_spherical(far_sample, query, "dn", "background", k = 1)
  expect_lt(near$dn, far$dn)
  out <- interpolate_spherical(samples, spherical_queries(90), "dn", "background")
  expect_equal(out$dn, 100)
})

test_that("rule is respected and zero filling weight disables the blend", {
  samples <- spherical_samples(c(0, 30), c(2, 20))
  query <- spherical_queries(5)
  any <- interpolate_spherical(samples, query, "dn", k = 2, rule = "any")
  all <- interpolate_spherical(samples, query, "dn", k = 2, rule = "all")
  expect_equal(any$dn, 2)
  expect_true(is.na(all$dn))
  all_filled <- interpolate_spherical(samples, query, "dn", "background",
                                       k = 2, rule = "all")
  expect_equal(all_filled$dn, 100)
  off <- interpolate_spherical(samples, query, "dn", "background", 0,
                                 k = 2, rule = "all")
  expect_true(is.na(off$dn))
  query <- spherical_queries(c(4, 90))
  query$background <- NA_real_
  off <- interpolate_spherical(samples, query, "dn", "background", 0,
                                 search_radius = NULL, min_dist = NULL)
  pure <- interpolate_spherical(samples, query, "dn", search_radius = NULL)
  expect_identical(off, pure)
  active <- interpolate_spherical(samples, query, "dn", "background")
  expect_true(all(is.na(active$dn)))
})

test_that("missing samples are omitted before searching and undefined queries stay NA", {
  samples <- spherical_samples(c(0, 1, 2), c(NA, 10, 30))
  query <- spherical_queries(c(0, NA))
  out <- interpolate_spherical(samples, query, "dn", k = 2, rule = "all")
  expect_equal(out$dn, c(14, NA))
  out <- interpolate_spherical(samples, query, "dn", "background")
  expect_true(is.finite(out$dn[1]))
  expect_true(is.na(out$dn[2]))
  query <- spherical_queries(0)
  query$azimuth <- NA_real_
  samples$azimuth[1] <- NA_real_
  out <- interpolate_spherical(samples, query, "dn", k = 2)
  expect_equal(out$dn, 14)
  empty <- samples[FALSE, ]
  expect_true(is.na(interpolate_spherical(empty, query, "dn")$dn))
  expect_equal(interpolate_spherical(empty, query, "dn", "background")$dn, 100)
  out <- interpolate_spherical(samples, query[FALSE, ], "dn")
  expect_equal(nrow(out), 0)
  expect_identical(out$dn, numeric())
})

test_that("spherical distances handle the azimuth seam and both frames", {
  samples <- data.frame(zenith = pi / 2, azimuth = c(355, 5) * pi / 180,
                        dn = c(10, 30))
  query <- data.frame(zenith = pi / 2, azimuth = 0)
  out <- interpolate_spherical(samples, query, "dn", k = 2)
  expect_equal(out$dn, 20)
  samples$azimuth <- (samples$azimuth + .72) %% (2 * pi)
  query$azimuth <- (query$azimuth + .72) %% (2 * pi)
  expect_equal(interpolate_spherical(samples, query, "dn", k = 2)$dn, out$dn)
  columns <- c("off_axis", "internal_azimuth")
  names(samples)[1:2] <- names(query) <- columns
  internal <- interpolate_spherical(samples, query, "dn", k = 2,
                                      angle_columns = columns)
  expect_equal(internal$dn, out$dn)
})

test_that("invalid contracts fail with explicit errors", {
  samples <- spherical_samples()
  query <- spherical_queries()
  run <- function(...) interpolate_spherical(samples, query, "dn", ...)
  expect_error(run(filling_source = "missing"), "columns")
  expect_error(run(filling_source = "background", min_dist = 20), "smaller")
  expect_error(run(filling_source = "background", search_radius = NULL), "NULL")
  expect_error(run(filling_weight = -1), "non-negative")
  expect_error(run(p = 0), "positive")
  expect_error(run(search_radius = 181), "180")
  expect_error(run(k = NULL, search_radius = NULL), "Supply")
  expect_error(run(rule = "other"))
  expect_error(run(angle_columns = c("zenith", "zenith")), "different")
  query$dn <- 1
  expect_error(run(), "reserved")
  query$dn <- NULL
  samples$dn[1] <- Inf
  expect_error(run(), "Inf")
  samples$dn[1] <- 1
  query$background[1] <- Inf
  expect_error(run(filling_source = "background"), "Inf")
})

test_that("parallel interpolation agrees with sequential evaluation", {
  skip_on_cran()
  skip_if_not_installed("doParallel")
  skip_if(is.na(parallel::detectCores()) || parallel::detectCores() < 2)
  on.exit(foreach::registerDoSEQ(), add = TRUE)
  samples <- spherical_samples()
  query <- spherical_queries()
  # Background blending exercises IDW and the additional query data together.
  expected <- interpolate_spherical(samples, query, "dn", filling_source = "background")
  actual <- interpolate_spherical(samples, query, "dn", filling_source = "background",
                                   parallel = TRUE, cores = 2)
  expect_equal(actual, expected)
})
