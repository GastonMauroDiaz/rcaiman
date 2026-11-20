spherical_test_d <- function(zenith, azimuth, nrow = 2) {
  d <- terra::rast(nrows = nrow, ncols = length(zenith) / nrow, nlyrs = 2)
  terra::values(d) <- cbind(zenith, azimuth)
  names(d) <- c("zenith", "azimuth")
  d
}

test_that("equal-area labels belong to the sky and survive partial coverage", {
  zenith <- rep(c(0, 15, 45, 75, 90), each = 12)
  azimuth <- rep(seq(0, 330, 30), 5)
  d <- spherical_test_d(zenith, azimuth)
  full <- beckers_segmentation(d, 20)
  cropped <- terra::crop(d, terra::ext(-108, 156, 0, 90))
  expect_equal(terra::values(beckers_segmentation(cropped, 20)),
               terra::values(terra::crop(full, terra::ext(cropped))))
  masked <- terra::deepcopy(d)
  masked[azimuth < 180] <- NA
  partial <- beckers_segmentation(masked, 20)
  expect_equal(terra::values(partial)[azimuth >= 180],
               terra::values(full)[azimuth >= 180])
  expect_true(all(terra::values(partial)[azimuth < 180] == 0))
  expect_equal(attr(partial, "n_cells"), 20)
})

test_that("equal-area boundaries and zenith cap do not require full azimuth coverage", {
  d <- spherical_test_d(c(0, 10, 80, 90, 91, NA), c(NA, NA, NA, 0, 20, 30))
  seg <- beckers_segmentation(d, 20)
  expect_equal(as.numeric(terra::values(seg))[1:3], c(1, 1, 0))
  expect_gt(terra::values(seg)[4], 0)
  expect_equal(as.numeric(terra::values(seg))[5:6], c(0, 0))
  expect_equal(as.numeric(terra::values(beckers_segmentation(d, 1))), c(1,1,1,1,0,0))
  expect_error(beckers_segmentation(d, 20, shape_control = 1e6), "decreasing")
  expect_error(beckers_segmentation(d, 20, min_k0 = 0), "positive")
})

test_that("HEALPix preserves IDs across coverage and masks below-horizon directions", {
  d <- spherical_test_d(c(0, 30, 50, 90, 100, NA), c(NA, 45, 225, 0, 20, 30))
  seg <- healpix_segmentation(d, 4)
  expect_equal(as.numeric(terra::values(seg))[c(1,5,6)], c(0,0,0))
  expect_true(all(terra::values(seg)[2:4] > 0))
  expect_equal(attr(seg, "cell_solid_angle"), 4*pi/(12*4^2))
  copy <- terra::deepcopy(d)
  copy[2] <- NA
  expect_equal(terra::values(healpix_segmentation(copy, 4))[3:4], terra::values(seg)[3:4])
  expect_error(healpix_segmentation(d, 3), "power of two")
  expect_error(healpix_segmentation(d, 16384), "8192")
})

test_that("adaptive and fixed HEALPix share the same angular reference", {
  grid <- expand.grid(zenith = seq(5, 85, 10), azimuth = seq(5, 355, 10))
  d <- spherical_test_d(grid$zenith, grid$azimuth)
  r <- d$zenith * 0 + 1
  fixed <- as.numeric(terra::values(healpix_segmentation(d, 1)))
  tree <- as.numeric(terra::values(adaptive_healpix_segmentation(r, d)))
  mapping <- unique(data.frame(fixed, tree))
  expect_false(anyDuplicated(mapping$fixed) > 0)
  expect_false(anyDuplicated(mapping$tree) > 0)
  # No missing radiometry may influence the variability of a region.
  r[1] <- NA
  tree <- adaptive_healpix_segmentation(r, d)
  expect_equal(terra::values(tree)[1], 0)
  expect_true(all(terra::values(tree)[-1] > 0))
  expect_error(adaptive_healpix_segmentation(r, d, max_splittings = 14), "13")
})

test_that("empty and below-horizon coverage yields an empty segmentation", {
  for (z in list(rep(NA_real_, 6), rep(100, 6))) {
    d <- spherical_test_d(z, rep(30, 6))
    r <- d[[1]] * 0 + 1
    expect_equal(as.numeric(terra::values(beckers_segmentation(d, 20))), rep(0,6))
    expect_equal(as.numeric(terra::values(healpix_segmentation(d))), rep(0,6))
    expect_equal(as.numeric(terra::values(adaptive_healpix_segmentation(r, d))), rep(0,6))
  }
})


test_that("adaptive subdivisions remain within their fixed HEALPix parents", {
  grid <- expand.grid(zenith = seq(1,89,2), azimuth = seq(1,359,2))
  d <- spherical_test_d(grid$zenith, grid$azimuth)
  r <- d$zenith * 0
  set.seed(48)
  terra::values(r) <- runif(terra::ncell(r))
  tree <- as.numeric(terra::values(adaptive_healpix_segmentation(r,d,
    max_splittings=3,min_size_px=2)))
  fixed <- as.numeric(terra::values(healpix_segmentation(d,1)))
  expect_gt(length(unique(tree)), length(unique(fixed)))
  mapping <- unique(data.frame(tree,fixed))
  expect_false(anyDuplicated(mapping$tree) > 0)
  r[1:10] <- NA
  masked <- terra::deepcopy(d)
  masked[1:10] <- NA
  expect_equal(terra::values(adaptive_healpix_segmentation(r,d)),
               terra::values(adaptive_healpix_segmentation(r,masked)))
})
