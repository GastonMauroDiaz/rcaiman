sky_test_d <- function(zenith, azimuth) {
  d <- terra::rast(nrows = 3, ncols = 4, nlyrs = 2,
    xmin = 10, xmax = 18, ymin = -2, ymax = 4)
  terra::values(d) <- cbind(zenith, azimuth)
  names(d) <- c("zenith", "azimuth")
  d
}

test_that("rings use external zenith, including boundaries and partial coverage", {
  d <- sky_test_d(c(0, 30, 30.1, 60, 90, 90.1, 180, NA, 15, 45, 75, NA),
                 rep(NA_real_, 12))
  rings <- ring_segmentation(d, 30)
  expect_equal(as.numeric(terra::values(rings)), c(1,1,2,2,3,0,0,0,1,2,3,0))
  expect_true(terra::compareGeom(rings, d[[1]]))
  angles <- ring_segmentation(d, 30, return = "angle")
  expect_equal(as.numeric(terra::values(angles)), c(15,15,45,45,75,0,0,0,15,45,75,0))
  expect_identical(attr(angles, "return"), "angle")
  d$azimuth <- Inf
  expect_equal(terra::values(ring_segmentation(d, 30)), terra::values(rings))
  expect_error(ring_segmentation(d, 7), "whole number")
})

test_that("sectors ignore zenith and accept widths that divide 360", {
  d <- sky_test_d(rep(NA_real_, 12),
                 c(0, 60, 60.1, 120, 180, 240, 300, 359.9, NA, 30, 90, 270))
  sectors <- sector_segmentation(d, 60)
  expect_equal(as.numeric(terra::values(sectors)), c(1,1,2,2,3,4,5,6,0,1,2,5))
  expect_true(terra::compareGeom(sectors, d[[1]]))
  d$zenith <- Inf
  expect_equal(terra::values(sector_segmentation(d, 60)), terra::values(sectors))
  d$azimuth <- 360
  expect_error(sector_segmentation(d, 60), "360")
})

test_that("sky grids require both angles except in an unsubdivided first ring", {
  d <- sky_test_d(c(0, 15, 45, 75, 90, 100, NA, 15, 45, 75, 15, 30),
                 c(NA, 15, 45, 75, 0, 10, 10, NA, NA, 359, 345, 30))
  g <- sky_grid_segmentation(d, 30)
  expect_equal(as.numeric(terra::values(g)),
    c(0,1001,2002,3003,1003,0,0,0,0,12003,12001,1001))
  joined <- sky_grid_segmentation(d, 30, first_ring_different = TRUE)
  expect_equal(as.numeric(terra::values(joined)),
    c(1000,1000,2002,3003,1003,0,0,1000,0,12003,1000,1000))
  expect_equal(attr(g, "angle_width"), 30)
  expect_false(attr(g, "first_ring_different"))
  expect_error(sky_grid_segmentation(d, 90 / 1000), "1000")
  expect_error(sky_grid_segmentation(d, 30, NA), "NA")
})

test_that("empty coverage and malformed direction stacks have clear behavior", {
  d <- sky_test_d(rep(NA_real_, 12), rep(NA_real_, 12))
  for (fun in list(ring_segmentation, sector_segmentation, sky_grid_segmentation)) {
    expect_equal(as.numeric(terra::values(fun(d, 30))), rep(0, 12))
    expect_error(fun(d[[1]], 30), "two layers")
  }
  d$zenith <- 181
  expect_error(ring_segmentation(d, 30), "180")
})

test_that("nominal centers map to unique matching sky-grid labels", {
  centers <- sky_grid_centers(30)
  d <- terra::rast(nrows = 3, ncols = 12, nlyrs = 2)
  terra::values(d) <- as.matrix(centers) * 180 / pi
  names(d) <- c("zenith", "azimuth")
  labels <- as.numeric(terra::values(sky_grid_segmentation(d, 30)))
  expect_equal(labels, rep(1:12, 3) * 1000 + rep(1:3, each = 12))
  expect_equal(length(unique(labels)), nrow(centers))
})
