sky_mask_fixture <- function() {
  d <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 3, nlyrs = 2))
  terra::values(d) <- cbind(c(0, 0, 75, 75, 100, 100), rep(90, 6))
  names(d) <- c("zenith", "azimuth")
  seg <- beckers_segmentation(d, 20)
  cv <- terra::rast(d[[1]])
  terra::values(cv) <- c(1, 3, 5, 9, NA, NA)
  names(cv) <- "arbitrary_cv_name"
  bin <- function(values) terra::setValues(terra::rast(cv), as.logical(values))
  list(seg = seg, w = terra::setValues(terra::rast(seg), c(1, 1, 3, 3, NA, NA)), cv = cv,
    indices = data.frame(col = c(1, 2, 3, 1), row = c(1, 1, 1, 2), tag = NA),
    masks = list(restrictive = bin(c(1, 1, 0, 0, 0, 0)),
                 permissive = bin(c(1, 1, 1, 1, 0, 0))), bin = bin)
}

test_that("mask optimization compares common samples and exposes raw metrics", {
  f <- sky_mask_fixture()
  result <- optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks,
                              coverage_weight = .75, return = "all")
  expect_identical(result$mask, f$masks[[2]])
  expect_equal(result$index, 2L)
  expect_equal(result$data$name, names(f$masks))
  expect_equal(result$data$n_samples, c(2, 4))
  expect_equal(result$data$consistency_raw, c(2, 4))
  expect_equal(result$data$coverage_raw, c(.25, 1))
  expect_equal(result$data$ccm, c(.75, .25))
  expect_identical(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks, .25), f$masks[[1]])
  expect_identical(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks, 0), f$masks[[1]])
  expect_identical(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks, 1), f$masks[[2]])
})

test_that("constant criteria do not prevent discrimination by another criterion", {
  f <- sky_mask_fixture()
  masks <- list(f$bin(c(1, 0, 1, 0, 0, 0)), f$bin(c(0, 1, 0, 1, 0, 0)))
  result <- optimize_sky_mask(f$cv, f$indices, f$seg, f$w, masks, return = "all")
  expect_equal(result$index, 1L)
  expect_equal(result$data$coverage, c(0, 0))
  expect_equal(result$data$consistency, c(0, 1))
  f$cv[] <- 1
  result <- optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks, return = "all")
  expect_equal(result$index, 2L)
  expect_equal(result$data$consistency, c(0, 0))
  expect_warning(result <- optimize_sky_mask(f$cv, f$indices, f$seg, f$w, masks,
    return = "all"), "same score")
  expect_equal(result$index, 1L)
  expect_identical(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, masks[1]), masks[[1]])
})

test_that("empty candidates are reported and excluded from normalization", {
  f <- sky_mask_fixture()
  masks <- c(list(empty = f$bin(rep(0, 6))), f$masks)
  expect_warning(result <- optimize_sky_mask(f$cv, f$indices, f$seg, f$w, masks,
    .75, return = "all"), "excluded: 1")
  expect_equal(result$index, 3L)
  expect_equal(result$data$n_samples[1], 0)
  expect_equal(result$data$coverage_raw[1], 0)
  expect_true(is.na(result$data$consistency_raw[1]))
  expect_true(is.na(result$data$ccm[1]))
  expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, masks[1]), "No candidate")
})

test_that("repeated indices affect the median, not coverage", {
  f <- sky_mask_fixture()
  repeated <- f$indices[c(1, 1, 1, 2, 3, 4), ]
  result <- optimize_sky_mask(f$cv, repeated, f$seg, f$w, f$masks, .75, return = "all")
  expect_equal(result$data$n_samples, c(4, 6))
  expect_equal(result$data$consistency_raw, c(1, 2))
  expect_equal(result$data$coverage_raw, c(.25, 1))
})

test_that("invalid selected data and incompatible grids are rejected", {
  f <- sky_mask_fixture()
  expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, list()), "nonempty")
  expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks, 1.1), "coverage_weight")
  expect_error(optimize_sky_mask(f$cv, f$indices[FALSE, ], f$seg, f$w, f$masks), "at least one row")
  expect_error(optimize_sky_mask(f$cv, data.frame(col = 3, row = 2), f$seg, f$w, f$masks),
               "positive cell IDs")
  for (value in c(NA, Inf, -1)) {
    f$cv[1] <- value
    expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks), "selected CV values")
  }
  f <- sky_mask_fixture()
  f$masks[[1]][1] <- NA
  expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks), "NA-free")
  f <- sky_mask_fixture()
  terra::ext(f$seg) <- terra::ext(0, 30, 0, 20)
  expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks), "same geometry")
})


test_that("mask candidates share the restricted reference area", {
  f <- sky_mask_fixture()
  m <- f$bin(c(1, 1, 1, 0, 0, 0))
  indices <- f$indices[1:3, ]
  result <- optimize_sky_mask(f$cv, indices, f$seg, f$w, f$masks,
    coverage_weight = .75, return = "all", m = m)
  expect_equal(result$data$coverage_raw, c(2/5, 1))
  for (i in seq_along(f$masks)) {
    selected <- rem_unselected(indices, f$masks[[i]])
    expected <- assess_sampling_uniformity(selected, f$seg, f$w, m)
    expect_equal(result$data$coverage_raw[i], expected$coverage)
  }
  expect_error(optimize_sky_mask(f$cv, f$indices, f$seg, f$w, f$masks, m = m),
               "reference domain")
})
