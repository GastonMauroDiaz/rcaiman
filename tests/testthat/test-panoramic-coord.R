test_that("panoramic coordinates retain padding values and row order", {
  sampling_points <- data.frame(
    zen = c(5, 30, 60, 85) * pi / 180,
    azi = c(5, 355, 150, 200) * pi / 180,
    dn = c(1, 2, 3, 4),
    x = 99,
    y_down = 99,
    col = 100,
    row = 100,
    y = 99,
    relative_radius = 0.5,
    image_angle = 0
  )
  observed <- to_panoramic(sampling_points, angle_width = 10, padding_deg = 10)
  sampling_coord <- observed$sampling_coord
  # Values recorded before migration; only the spatial names are different.
  expect_named(observed, c("sampling_coord", "raster"))
  expect_named(sampling_coord, c("x", "y_down", "zen", "azi", "dn"))
  expect_equal(sampling_coord$x, c(1.5, 36.5, 16, 21, 1.5, 37.5, 37.5, 0.5))
  expect_equal(sampling_coord$y_down, c(1.5, 4, 7, 9.5, 0.5, 1.5, 0.5, 4))
  source_rows <- c(1, 2, 3, 4, 1, 1, 1, 2)
  expect_equal(sampling_coord$dn, sampling_points$dn[source_rows])
  expect_equal(sampling_coord$zen, sampling_points$zen[source_rows])
  expect_equal(sampling_coord$azi, sampling_points$azi[source_rows])
  expect_equal(dim(observed$raster), c(10, 38, 1))
  expect_equal(unname(as.vector(terra::ext(observed$raster))), c(0, 38, 0, 10))
})

test_that("unpadded panoramic positions remain continuous including edges", {
  sampling_points <- data.frame(
    zen = c(0, 22.5, 90, 22.5) * pi / 180,
    azi = c(0, 12.5, 359.9, 12.5) * pi / 180,
    id = c("zenith", "inside", "horizon", "duplicate")
  )
  observed <- to_panoramic(sampling_points, angle_width = 5, padding_deg = 0)
  expect_equal(observed$sampling_coord$x, c(0, 2.5, 71.98, 2.5))
  expect_equal(observed$sampling_coord$y_down, c(0, 4.5, 18, 4.5))
  expect_identical(observed$sampling_coord$id, sampling_points$id)
  empty <- to_panoramic(sampling_points[FALSE, ], 5, 0)
  expect_equal(nrow(empty$sampling_coord), 0)
  expect_named(empty$sampling_coord, c("x", "y_down", "zen", "azi", "id"))
})

test_that("panoramic return selects cells containing continuous coordinates", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  angles <- data.frame(zen = 0, azi = 0)
  for (angle_width in c(5, 10)) {
    pano <- to_panoramic(angles, angle_width, padding_deg = 0)$raster
    terra::values(pano) <- seq_len(terra::ncell(pano))
    observed <- from_panoramic(pano, z, a, angle_width)
    valid <- which(!is.na(terra::values(z)))
    # Independent terra lookup checks both the half-cell shift and y flip.
    target_xy <- cbind(
      x = terra::values(a)[valid] / angle_width,
      y = nrow(pano) - terra::values(z)[valid] / angle_width
    )
    expected_cells <- terra::cellFromXY(pano, target_xy)
    expect_equal(terra::values(observed)[valid], as.numeric(expected_cells))
    expect_true(all(is.na(terra::values(observed)[-valid])))
    expect_true(terra::compareGeom(observed, z))
  }
})

test_that("panoramic coord composes with planar interpolation and padding", {
  z <- zenith_image(10, equidistant_coef)
  a <- azimuth_image(z)
  sampling_indices <- raster_grid_in_circle(z, a, spacing = 2)
  sampling_indices$dn <- 7
  pano <- to_panoramic(sampling_indices, angle_width = 5, padding_deg = 10)
  interpolated <- interpolate_planar(
    pano$sampling_coord, r = pano$raster,
    search_radius = 100, parallel = FALSE
  )
  unpadded <- rem_panoramic_padding(interpolated, 5, 10)
  observed <- from_panoramic(unpadded, z, a, 5)
  expect_equal(dim(unpadded), c(18, 72, 1))
  expect_equal(unname(as.vector(terra::ext(unpadded))), c(0, 72, 0, 18))
  expect_equal(terra::values(observed)[!is.na(terra::values(z))],
               rep(7, sum(!is.na(terra::values(z)))))
})
