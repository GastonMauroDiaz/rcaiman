test_that("raster_grid_in_circle() requires integer-like spacing", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)

  expect_error(
    raster_grid_in_circle(z, a, spacing = 1.5),
    "spacing.*integer-like"
  )

  observed <- raster_grid_in_circle(z, a, spacing = 2)
  expect_true(all(observed$col == round(observed$col)))
  expect_true(all(observed$row == round(observed$row)))
})
