test_that("sky-grid centers use external angular midpoints and stable row order", {
  grid <- sky_grid_centers(30)
  expect_identical(names(grid), c("zenith", "azimuth"))
  expect_equal(nrow(grid), 36L)
  expect_equal(grid$zenith, rep(c(15, 45, 75), each = 12) * pi / 180)
  expect_equal(grid$azimuth, rep(seq(15, 345, 30), 3) * pi / 180)
  expect_identical(grid, sky_grid_centers(30))
  expect_null(attr(grid, "out.attrs"))
  fractional <- sky_grid_centers(90 / 7)
  expect_equal(nrow(fractional), 7 * 28)
  expect_true(all(fractional$zenith > 0 & fractional$zenith < pi / 2))
  expect_true(all(fractional$azimuth > 0 & fractional$azimuth < 2 * pi))
  expect_error(sky_grid_centers(7), "whole number")
  expect_error(sky_grid_centers(120), "whole number")
  expect_error(sky_grid_centers(0), "positive")
  expect_error(sky_grid_centers(NA_real_), "NA")
  expect_error(sky_grid_centers(Inf), "Inf")
})

test_that("sky-grid queries locate the same sky cells in differently oriented images", {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(200, 150), normalization_radius = 140,
    normalization_angle = 90)
  g <- set_internal_geometry(c(400, 300), geometry)
  queries <- sky_grid_centers(15)
  for (elevation in c(90, 32.5)) {
    d <- direction_image(g, elevation = elevation)
    coord <- reproject_points(queries, geometry, elevation = elevation)
    indices <- add_indices(coord)
    inside <- with(indices, !is.na(col) & !is.na(row) &
      col >= 1 & col <= ncol(d) & row >= 1 & row <= nrow(d))
    indices <- rem_unselected(indices[inside, ], !is.na(d$zenith))
    expect_gt(nrow(indices), 0)
    seg <- sky_grid_segmentation(d, 15)
    labels <- extract_dn(seg, indices, use_window = FALSE)
    expected <- ceiling(indices$azimuth * 180 / pi / 15) * 1000 +
      ceiling(indices$zenith * 180 / pi / 15)
    expect_equal(labels[[ncol(labels)]], expected)
    measured <- add_directions(indices[c("col", "row")], d = d)
    expect_equal(measured$zenith, indices$zenith, tolerance = .02)
  }
})
