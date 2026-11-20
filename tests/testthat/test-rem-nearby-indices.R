test_that("rem_nearby_indices() filters in cell-index space", {
  sampling_indices <- data.frame(
    col = c(1, 2, 10),
    row = c(1, 1, 10),
    id = letters[1:3]
  )

  observed <- rem_nearby_indices(sampling_indices, min_dist = 2)

  expect_equal(observed, sampling_indices[c(1, 3), ])
})


test_that("rem_nearby_indices() filters in angular space", {
  sampling_indices <- data.frame(
    col = 1:3,
    row = 1:3,
    zenith = rep(pi / 4, 3),
    azimuth = c(0, pi / 180, pi)
  )

  observed <- rem_nearby_indices(
    sampling_indices,
    min_dist = 10,
    space = "spherical"
  )

  expect_equal(observed, sampling_indices[c(1, 3), ])
})

test_that("ranking controls retention while output preserves input order", {
  x <- data.frame(col = c(1, 2, 10), row = c(1, 1, 10),
                  dn = c(100, 200, 300), tag = letters[1:3],
                  zenith = rep(pi / 4, 3), azimuth = c(0, pi / 180, pi))
  for (space in c("planar", "spherical")) {
    expect_equal(rem_nearby_indices(x, 2, "dn", "higher", space), x[c(2, 3), ])
    expect_equal(rem_nearby_indices(x, 2, 3, "lower", space), x[c(1, 3), ])
    x$dn[2] <- 100
    expect_equal(rem_nearby_indices(x, 2, "dn", "higher", space), x[c(1, 3), ])
    x$dn[2] <- 200
  }
})

test_that("angular selection agrees across rotated internal and external frames", {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(20, 15), normalization_radius = 15,
    normalization_angle = 90)
  g <- set_internal_geometry(c(40, 30), geometry)
  d <- direction_image(g, sensor_rotation = 17, azimuth_rotation = 42,
                       elevation = 32.5)
  indices <- expand.grid(col = seq(1, 40, 3), row = seq(1, 30, 3))
  indices <- rem_unselected(indices, !is.na(g$off_axis))
  x <- add_directions(indices, g = g, d = d)
  internal <- rem_nearby_indices(x, 17, space = "spherical",
    angle_columns = c("off_axis", "internal_azimuth"))
  external <- rem_nearby_indices(x, 17, space = "spherical")
  expect_equal(external, internal)
  xyz <- as.matrix(.spherical_neighbor_vectors(external, c("zenith", "azimuth")))
  angle <- acos(pmax(pmin(tcrossprod(xyz), 1), -1))
  expect_true(all(angle[lower.tri(angle)] >= 17 * pi / 180 - 1e-12))
})

test_that("zero distance and empty inputs preserve rows and columns", {
  x <- data.frame(col = c(1, 1), row = c(1, 1), dn = c(1, 2),
                  zenith = c(.5, .5), azimuth = c(1, 1))
  for (space in c("planar", "spherical")) {
    expect_equal(rem_nearby_indices(x, 0, "dn", "higher", space), x)
    expect_equal(rem_nearby_indices(x, 1, space = space), x[1, ])
    expect_equal(rem_nearby_indices(x[FALSE, ], 1, space = space), x[FALSE, ])
    expect_equal(rem_nearby_indices(x[1, ], 1, space = space), x[1, ])
  }
})

test_that("polar directions allow missing azimuth but other missing directions fail", {
  x <- data.frame(col = 1:3, row = 1:3, zenith = c(0, .01, pi),
                  azimuth = c(NA, .5, NA))
  expect_equal(rem_nearby_indices(x, 2, space = "spherical"), x[c(1, 3), ])
  x$zenith[1] <- .5
  expect_error(rem_nearby_indices(x, 2, space = "spherical"), "must not be missing")
  x$zenith[1] <- NA_real_
  expect_error(rem_nearby_indices(x, 2, space = "spherical"), "must not be missing")
})

test_that("invalid indices, rankings and angular parameters are rejected", {
  x <- data.frame(col = 1:3, row = 1:3, dn = 1:3,
                  zenith = rep(.5, 3), azimuth = c(0, 1, 2))
  expect_error(rem_nearby_indices(x, 181, space = "spherical"), "180")
  expect_error(rem_nearby_indices(x, space = "spherical",
    angle_columns = c("zenith", "zenith")), "different columns")
  x$dn[1] <- NA_real_
  expect_error(rem_nearby_indices(x, col_id = "dn", prioritize = "higher"), "ranking values")
  x$dn[1] <- Inf
  expect_error(rem_nearby_indices(x, col_id = "dn", prioritize = "lower"), "ranking values")
  x$col[1] <- 1.5
  expect_error(rem_nearby_indices(x), "sampling_indices\\$col")
  x$col[1] <- 0
  expect_error(rem_nearby_indices(x), "sampling_indices\\$col")
})
