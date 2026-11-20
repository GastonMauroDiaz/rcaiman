# Protect against incidental algorithmic changes within the intended use of
# the package. References were recorded before the coord migration on
# 2026-09-04, then supplied through coord_from_indices() at the same centers.
# This rectangular raster is an ordinary image template, not a fabricated
# zenith or azimuth image. Angular fixtures must use the package builders.

planar_coord_reference <- function() {
  r <- set_rcaiman_geometry(terra::rast(nrows = 5, ncols = 7))
  terra::values(r) <- 0
  names(r) <- "dn"
  sampling_coord <- coord_from_indices(
    col = c(1, 7, 1, 7, 2, 5, 3, 6),
    row = c(1, 1, 5, 5, 3, 2, 4, 4)
  )
  sampling_coord$dn <- c(2, 7, 4, 11, 3, 9, 6, 10)
  list(r = r, sampling_coord = sampling_coord)
}

test_that("both IDW engines retain their cell-center reference", {
  fixture <- planar_coord_reference()
  expected <- c(
    2, 2.703703703704, 4.461538461538, 6.924528301887, 8.703703703704,
    7.827586206897, 7, 2.615384615385, 3.058823529412, 5.25,
    7.714285714286, 9, 8.571428571429, 7.758620689655, 3, 3,
    4.909090909091, 6.6, 8.941176470588, 9.352941176471, 9.5,
    4, 4.4, 6, 6.275862068966, 9.166666666667, 10, 10.411764705882,
    4, 4.428571428571, 5.241379310345, 6.515151515152, 9.421052631579,
    10.428571428571, 11
  )
  for (engine in c("cpp", "terra")) {
    observed <- interpolate_planar(
      fixture$sampling_coord, r = fixture$r,
      k = 3, p = 2, search_radius = 20, engine = engine, parallel = FALSE
    )
    expect_equal(as.vector(terra::values(observed)), expected,
                 tolerance = 1e-10)
    expect_true(terra::compareGeom(observed, fixture$r))
  }
})

test_that("IDW uses subpixel positions without collapsing a shared cell", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 2))
  sampling_coord <- data.frame(
    x = c(0.25, 1.25), y_down = c(0.5, 0.5), dn = c(0, 100)
  )
  for (engine in c("cpp", "terra")) {
    observed <- interpolate_planar(
      sampling_coord, r = r, k = 2, p = 2, search_radius = 3,
      engine = engine, parallel = FALSE
    )
    expect_equal(as.vector(terra::values(observed)), c(10, 2500 / 26))

    same_cell_coord <- sampling_coord
    same_cell_coord$x <- c(0.25, 0.75)
    observed <- interpolate_planar(
      same_cell_coord, r = r, k = 2, search_radius = 3,
      engine = engine, parallel = FALSE
    )
    expect_equal(unname(terra::values(observed)[1, 1]), 50)
  }
})

test_that("IDW omits missing samples and returns NA for an empty sample set", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 3))
  points <- data.frame(x = c(0.5, 1.5, NA, Inf),
                       y_down = c(0.5, 0.5, 1, 1), dn = c(7, NA, 20, 30))
  for (engine in c("cpp", "terra")) {
    actual <- interpolate_planar(points, r = r, engine = engine,
                                  search_radius = 10, parallel = FALSE)
    expect_equal(as.vector(terra::values(actual)), rep(7, 6))
    empty <- interpolate_planar(points[FALSE, ], r = r, engine = engine,
                                 parallel = FALSE)
    expect_true(all(is.na(terra::values(empty))))
    expect_true(terra::compareGeom(empty, r))
  }
})

test_that("triangulation retains the rasterized-plane reference", {
  fixture <- planar_coord_reference()
  observed <- triangulate_planar(fixture$sampling_coord, r = fixture$r)
  # The missing leftmost column is an existing rasterization edge effect.
  # This reference records it; the coord migration does not repair it.
  expected <- c(
    NaN, 2.833333581686, 3.666666895151, 4.500000208616, 5.333333522081,
    6.166666835546, 7.000000149012, NaN, 3.428571358323, 5.285714164376,
    7.142856970429, 8.999999776483, 8.399999827147, 8, NaN, 3, 5.25, 7.5,
    8.833333604038, 9.199999839067, 9, NaN, 4.333333134651, 5.999999761581,
    7.333333559334, 8.666666932404, 10.000000305474, 10, NaN, 5.166666522622,
    6.333333149552, 7.499999776483, 8.666666403413, 9.833333030343,
    10.999999657273
  )
  expect_equal(as.vector(terra::values(observed)), expected, tolerance = 1e-6)
  expect_true(terra::compareGeom(observed, fixture$r))

  terra::values(fixture$r) <- NA_real_
  empty_template <- triangulate_planar(
    fixture$sampling_coord, col_id = 3, r = fixture$r
  )
  expect_equal(terra::values(empty_template), terra::values(observed))
})

test_that("triangulation keeps continuous vertex positions", {
  fixture <- planar_coord_reference()
  sampling_coord <- fixture$sampling_coord
  sampling_coord$x <- sampling_coord$x + 0.1
  sampling_coord$y_down <- sampling_coord$y_down + 0.2
  observed <- triangulate_planar(
    sampling_coord, r = fixture$r, return = "all"
  )
  vertex_xy <- unique(terra::geom(observed$mesh)[, c("x", "y")])
  expected_xy <- cbind(
    x = sampling_coord$x,
    y = y_from_y_down(sampling_coord$y_down, nrow(fixture$r))
  )
  ordered_xy <- function(x) unname(x[order(x[, 1], x[, 2]), , drop = FALSE])
  expect_equal(ordered_xy(vertex_xy), ordered_xy(expected_xy))
  expect_true(any(is.finite(terra::values(observed$raster))))
})

test_that("planar interfaces do not silently reinterpret indices as coord", {
  fixture <- planar_coord_reference()
  sampling_indices <- data.frame(col = 1:3, row = 1:3, dn = 1:3)
  for (interpolator in list(interpolate_planar, triangulate_planar)) {
    expect_error(interpolator(sampling_indices, r = fixture$r), "y_down")
  }
  expect_error(fit_trend_surface(sampling_indices), "y_down")
})

test_that("modeled directions compose with reprojection and interpolation", {
  source_coord <- coord_from_indices(
    col = c(1, 4, 8, 2, 6, 7, 3, 5), row = c(1, 1, 2, 3, 4, 5, 6, 6)
  )
  source_coord$dn <- rep(7, nrow(source_coord))
  parameters <- list(
    off_axis_coef = 1,
    sin_internal_azimuth_coef = c(0, 1),
    cos_internal_azimuth_coef = c(-1, 0),
    atan_scale = 0.6
  )
  source_coord <- internal_directions_from_coord(source_coord,
    list(model = "anisotropic_projection", center_of_vision = c(4, 3),
         parameters = parameters), dim = c(8, 6))
  target_coord <- coord_from_internal_directions(source_coord, list(
    model = "radial_projection", center_of_vision = c(5, 5),
    parameters = equidistant_coef, normalization_radius = 5,
    normalization_angle = 90
  ))
  z <- set_rcaiman_geometry(terra::rast(nrows = 10, ncols = 10))

  expect_identical(target_coord$dn, source_coord$dn)
  for (engine in c("cpp", "terra")) {
    observed <- interpolate_planar(
      target_coord, r = z, search_radius = 20, engine = engine,
      parallel = FALSE
    )
    expect_equal(as.vector(terra::values(observed)), rep(7, terra::ncell(z)))
  }
})
