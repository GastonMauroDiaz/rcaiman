radial_standard <- function(center = c(4, 3), radius = 10) {
  list(model = "radial_projection", center_of_vision = center,
       parameters = 2 / pi, normalization_radius = radius,
       normalization_angle = 90)
}

test_that("radial point reprojection inverts a synthetic destination", {
  geometry <- list(model = "radial_projection", center_of_vision = c(3, 2),
                   parameters = c(0, 1 / (pi / 3)^2),
                   normalization_radius = 3.5, normalization_angle = 60,
                   max_off_axis_angle = 140)
  g <- set_internal_geometry(c(8, 6), geometry)
  d <- direction_image(g)
  indices <- expand.grid(col = seq_len(ncol(g)), row = seq_len(nrow(g)))
  indices <- rem_unselected(indices, !is.na(d[["zenith"]]))
  points <- add_directions(indices, d = d)
  points$id <- seq_len(nrow(points))
  points$value <- seq_len(nrow(points)) / 3
  before <- points
  actual <- reproject_points(points, geometry)
  valid <- !is.na(points$zenith)
  expected <- coord_from_indices(indices$col, indices$row)
  expect_equal(actual$x[valid], expected$x[valid], tolerance = 1e-8)
  expect_equal(actual$y_down[valid], expected$y_down[valid], tolerance = 1e-8)
  expect_identical(actual[c("zenith", "azimuth", "id", "value")],
                   points[c("zenith", "azimuth", "id", "value")])
  expect_identical(points, before)

  class(geometry) <- "geometry_spec"
  geometry$notes <- "A synthetic target, not an acquisition calibration."
  expect_identical(reproject_points(points, geometry), actual)
})

test_that("point reprojection retains missing and out-of-domain rows", {
  points <- data.frame(zenith = c(0, 45, 100, NA, 45, 0) * pi / 180,
                       azimuth = c(NA, NA, 90, 0, 270, 270) * pi / 180,
                       value = c(10, 20, 30, 40, 50, 60))
  result <- reproject_points(points, radial_standard())
  expect_equal(result$x, c(4, NA, NA, NA, 9, 4))
  expect_equal(result$y_down, c(3, NA, NA, NA, 3, 3))
  expect_identical(result[names(points)], points)
  empty <- points[FALSE, ]
  expect_equal(nrow(reproject_points(empty, radial_standard())), 0)
})

test_that("reprojection preserves aligned bands with a small search radius", {
  geometry <- radial_standard()
  g <- set_internal_geometry(c(8, 6), geometry)
  d <- direction_image(g)
  r <- g
  original <- cbind(0:47, 100:147)
  original[15, 1] <- NA_real_
  terra::values(r) <- original
  names(r) <- c("red", "nir")
  before <- terra::values(r)
  mask <- g[[1]] > 0
  mask[10] <- FALSE
  for (engine in c("cpp", "terra")) {
    output <- reproject_image(r, d, c(8, 6), geometry, m = mask,
                              search_radius = 0.1, engine = engine,
                              parallel = FALSE)
    expected <- original
    expected[10, ] <- NA_real_
    expect_equal(unname(terra::values(output)), expected)
    expect_identical(names(output), c("red", "nir"))
    expect_true(terra::compareGeom(r, output))
    expect_equal(unname(terra::values(output)[1, 1]), 0)
  }
  expect_identical(terra::values(r), before)
})

test_that("inclined partial images occupy their part of a larger standard", {
  targets <- radial_standard(c(20, 20), 20)
  sources <- list(
    radial_standard(c(2, 1), 30),
    list(model = "anisotropic_projection", center_of_vision = c(2, 1),
         parameters = list(off_axis_coef = 0.1,
                           sin_internal_azimuth_coef = c(0, 1),
                           cos_internal_azimuth_coef = c(-1, 0),
                           atan_scale = 0.6))
  )
  for (source in sources) {
    g <- set_internal_geometry(c(4, 2), source)
    d <- direction_image(g, elevation = 32.5)
    r <- g[[1]] * 0 + 7
    names(r) <- "radiance"
    output <- reproject_image(r, d, c(40, 40), targets,
                              search_radius = 1.5, parallel = FALSE)
    values <- terra::values(output, mat = FALSE)
    expect_equal(dim(output), c(40, 40, 1))
    expect_true(any(is.finite(values)))
    expect_equal(values[is.finite(values)], rep(7, sum(is.finite(values))),
                 tolerance = 1e-12)
    expect_true(all(is.na(values[c(1, 40, 781, 820, 1600)])))
    xy <- terra::xyFromCell(output, which(is.finite(values)))
    expect_true(all(40 - xy[, 2] > 28))
    expect_identical(names(output), "radiance")
  }
})

test_that("image reprojection forwards the planar interpolation controls", {
  geometry <- radial_standard(c(2, 2), 10)
  g <- set_internal_geometry(c(2, 1), radial_standard(c(1, 1), 10))
  d <- direction_image(g)
  # Place two observations at x = 0.25 and 1.25, y_down = 2.
  directions <- terra::values(d)
  directions[, "zenith"] <- c(1.75, 0.75) * 9
  directions[, "azimuth"] <- 90
  terra::values(d) <- directions
  r <- g[[1]]
  terra::values(r) <- c(0, 100)
  projected <- data.frame(x = c(0.25, 1.25), y_down = 2, dn = c(0, 100))
  template <- set_internal_geometry(c(4, 4), geometry)[[1]]
  # Exercise both controls together; their mathematics is tested in planar-coord.
  for (k in c(1, 2)) {
    p <- if (k == 1) 1 else 3
    expected <- interpolate_planar(projected, r = template, k = k, p = p,
                                    search_radius = 4, parallel = FALSE)
    actual <- reproject_image(r, d, c(4, 4), geometry, k = k, p = p,
                               search_radius = 4, parallel = FALSE)
    expect_equal(terra::values(actual), terra::values(expected),
                 ignore_attr = TRUE, tolerance = 1e-12)
  }
  insufficient <- reproject_image(r, d, c(4, 4), geometry, k = 3,
                                    rigorous = TRUE, parallel = FALSE)
  expect_true(all(is.na(terra::values(insufficient))))
})

test_that("a standard limits the target domain and handles absent samples", {
  geometry <- radial_standard(c(4, 3), 2)
  g <- set_internal_geometry(c(8, 6), geometry)
  d <- direction_image(g)
  r <- g[[1]] * 0 + 2
  out <- reproject_image(r, d, c(8, 6), geometry,
                         search_radius = 100, parallel = FALSE)
  outside <- is.na(terra::values(g[[1]], mat = FALSE))
  expect_true(all(is.na(terra::values(out, mat = FALSE)[outside])))
  expect_true(all(is.finite(terra::values(out, mat = FALSE)[!outside])))
  terra::values(r) <- NA_real_
  for (engine in c("cpp", "terra")) {
    out <- reproject_image(r, d, c(5, 7), geometry,
                           engine = engine, parallel = FALSE)
    expect_equal(dim(out), c(7, 5, 1))
    expect_true(all(is.na(terra::values(out))))
  }
})

test_that("reprojection validates angular units, grids and destination model", {
  geometry <- radial_standard()
  g <- set_internal_geometry(c(8, 6), geometry)
  d <- direction_image(g)
  points <- data.frame(zenith = pi / 6, azimuth = pi / 2)
  expect_error(reproject_points(data.frame(zen = 1, azi = 1), geometry),
                "zenith, azimuth", fixed = TRUE)
  bad <- points
  bad$zenith <- pi + 0.01
  expect_error(reproject_points(bad, geometry), "zenith")
  bad <- points
  bad$azimuth <- 2 * pi
  expect_error(reproject_points(bad, geometry), "azimuth")
  expect_error(reproject_image(g[[1]], g, c(8, 6), geometry), "zenith")
  cropped <- crop_caim(d, c(1, 1), 4, 3)
  expect_error(reproject_image(g[[1]], cropped, c(8, 6), geometry),
                "same geometry")
  anisotropic <- list(model = "anisotropic_projection", center_of_vision = c(4, 3),
                      parameters = list(off_axis_coef = 1, atan_scale = 0.5,
                                        sin_internal_azimuth_coef = c(0, 1),
                                        cos_internal_azimuth_coef = c(-1, 0)))
  expect_error(reproject_points(points, anisotropic), "radial_projection")
})

test_that("oriented point destinations invert direction_image continuously", {
  geometry <- radial_standard(c(3, 2), 4)
  geometry$max_off_axis_angle <- 150
  g <- set_internal_geometry(c(8, 6), geometry)
  indices <- expand.grid(col = seq_len(ncol(g)), row = seq_len(nrow(g)))
  orientations <- list(c(0, 0, 32.5), c(90, 90, 32.5), c(270, 180, 0),
                       c(17, 231, -24), c(-90, 450, 90), c(37, -120, -90))
  for (orientation in orientations) {
    d <- direction_image(g, orientation[1], orientation[2], orientation[3])
    retained <- rem_unselected(indices, !is.na(d[["zenith"]]))
    points <- add_directions(retained, g = g, d = d)
    expected <- coord_from_indices(retained$col, retained$row)
    out <- reproject_points(points, geometry,
                             orientation[1], orientation[2], orientation[3])
    valid <- !is.na(points$zenith)
    expect_equal(out$x[valid], expected$x[valid], tolerance = 1e-8)
    expect_equal(out$y_down[valid], expected$y_down[valid], tolerance = 1e-8)
    angles <- c("off_axis", "internal_azimuth", "zenith", "azimuth")
    expect_identical(out[angles], points[angles])
  }
})

test_that("destination rotations place known rays and use an off-axis limit", {
  geometry <- radial_standard(c(20, 20), 9)
  points <- data.frame(zenith = c(pi / 2, pi / 2, 0, pi),
                       azimuth = c(pi, 3 * pi / 2, NA, NA))
  horizontal <- reproject_points(points, geometry, elevation = 0)
  expect_equal(horizontal$x, c(20, 29, 20, 20))
  expect_equal(horizontal$y_down, c(20, 20, 11, 29))
  portrait <- reproject_points(points, geometry, sensor_rotation = 90,
                                elevation = 0)
  expect_equal(portrait$x, c(20, 20, 11, 29))
  expect_equal(portrait$y_down, c(20, 11, 20, 20))
  west <- reproject_points(points, geometry, azimuth_rotation = 90,
                            elevation = 0)
  expect_equal(west$x, c(11, 20, 20, 20))
  expect_equal(west$y_down, c(20, 20, 11, 29))

  geometry$max_off_axis_angle <- 70
  hinge <- data.frame(zenith = c(57.5, 0, 90, 120, 60, NA, 45) * pi / 180,
                      azimuth = c(pi, NA, pi, pi, 0, pi, NA))
  out <- reproject_points(hinge, geometry, elevation = 32.5)
  expect_equal(out$x, c(rep(20, 4), NA, NA, NA))
  expect_equal(out$y_down, c(20, 14.25, 23.25, 26.25, NA, NA, NA))

  geometry$max_off_axis_angle <- 180
  antipode <- reproject_points(data.frame(zenith = 122.5 * pi / 180, azimuth = 0),
                                geometry, elevation = 32.5)
  expect_true(is.na(antipode$x) && is.na(antipode$y_down))
  nadir <- reproject_points(points[4, ], geometry, elevation = -90)
  expect_equal(nadir$x, 20)
  expect_equal(nadir$y_down, 20)
})

test_that("image reprojection uses all destination orientation arguments", {
  geometry <- radial_standard()
  g <- set_internal_geometry(c(8, 6), geometry)
  d <- direction_image(g, sensor_rotation = 90, azimuth_rotation = 73,
                        elevation = 32.5)
  r <- g[[1]]
  terra::values(r) <- seq_len(terra::ncell(r))
  for (engine in c("cpp", "terra")) {
    out <- reproject_image(r, d, c(8, 6), geometry, sensor_rotation = 90,
                             azimuth_rotation = 73, elevation = 32.5,
                             search_radius = 0.1, engine = engine,
                             parallel = FALSE)
    expect_equal(terra::values(out), terra::values(r))
    # Turning only the destination sensor by 180 degrees turns the image.
    turned <- reproject_image(r, d, c(8, 6), geometry, sensor_rotation = 270,
                                azimuth_rotation = 73, elevation = 32.5,
                                search_radius = 0.1, engine = engine,
                                parallel = FALSE)
    expect_equal(terra::values(turned, mat = FALSE),
                   rev(terra::values(r, mat = FALSE)))
  }
})

test_that("destination orientation has the same validation as direction_image", {
  geometry <- radial_standard()
  points <- data.frame(zenith = pi / 6, azimuth = pi / 2)
  expect_error(reproject_points(points, geometry, elevation = 91), "elevation")
  g <- set_internal_geometry(c(8, 6), geometry)
  d <- direction_image(g)
  expect_error(reproject_image(g[[1]], d, c(8, 6), geometry, elevation = 91),
                "elevation")
})
