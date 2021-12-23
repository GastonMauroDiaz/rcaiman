test_that("ring_segmentation() and sector_segmentation() works", {
  z <- zenith_image(1490, lens())
  a <- azimuth_image(z)
  expect_equal(max(rings_segmentation(z, 15)[], na.rm = TRUE), 90/15)
  expect_equal(max(sectors_segmentation(a, 15)[], na.rm = TRUE), 360/15)
})

test_that("sky_grid_segmentation() works", {
  z <- zenith_image(1000, lens())
  a <- azimuth_image(z)
  expect_equal(max(sky_grid_segmentation(z, a, 10)[], na.rm = TRUE),
               36009)
})
