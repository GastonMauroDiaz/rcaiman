test_that("radius searches use directions from the scripting workflow", {
  geometries <- list(
    list(model = "radial_projection", parameters = equidistant_coef,
         center_of_vision = c(10, 8), normalization_radius = 11,
         normalization_angle = 90),
    list(model = "anisotropic_projection", center_of_vision = c(10, 8),
         parameters = list(off_axis_coef = c(0.9, 0.1), atan_scale = 0.6,
                           sin_internal_azimuth_coef = c(0, 1),
                           cos_internal_azimuth_coef = c(-1, 0)))
  )
  for (geometry in geometries) {
    g <- set_internal_geometry(c(20, 16), geometry)
    d <- direction_image(g, sensor_rotation = 90, azimuth_rotation = 123,
                         elevation = 32.5)
    indices <- expand.grid(col = 1:20, row = 1:16)
    indices <- rem_unselected(indices, !is.na(g$off_axis) & !is.na(g$internal_azimuth))
    points <- add_directions(indices, g, d)
    queries <- add_directions(data.frame(col = c(10, 14, 5), row = c(8, 8, 5)),
                              g, d)
    internal <- nn_spherical(sampling_directions = points,
                              query_directions = queries, k = NULL, search_radius = 18,
                              angle_columns = c("off_axis", "internal_azimuth"),
                              parallel = FALSE)
    external <- nn_spherical(points, queries, k = NULL, search_radius = 18,
                              parallel = FALSE)
    for (i in seq_len(nrow(queries))) {
      distance <- calc_spherical_distance(
        points$off_axis, points$internal_azimuth,
        queries$off_axis[i], queries$internal_azimuth[i])
      expected <- which(distance <= 18 * pi / 180)
      expect_setequal(internal[[i]]$index, expected)
      expect_equal(internal[[i]]$distance, distance[internal[[i]]$index],
                   tolerance = 1e-7)
      expect_true(all(diff(internal[[i]]$distance) >= 0))
      expect_setequal(external[[i]]$index, expected)
    }
    parallel <- nn_spherical(points, queries, k = NULL, search_radius = 18,
                              cores = 2, leave_free = 0)
    expect_identical(parallel, external)
  }
})

test_that("nearest, capped and all-within-radius searches have distinct contracts", {
  points <- data.frame(zenith = c(0, 10, 20, 30, 40) * pi / 180, azimuth = 0)
  query <- points[1, ]
  all <- nn_spherical(points, query, k = NULL, search_radius = 20,
                       parallel = FALSE)[[1]]
  expect_identical(all$index, 1:3)
  expect_equal(all$distance, c(0, 10, 20) * pi / 180)
  expect_identical(nn_spherical(points, query, k = 99, search_radius = 20,
                                parallel = FALSE)[[1]], all)
  expect_identical(nn_spherical(points, query, k = 2, search_radius = 20,
                                parallel = FALSE)[[1]]$index, 1:2)
  expect_length(nn_spherical(points, query, k = 4, search_radius = 20,
                             rule = "all", parallel = FALSE)[[1]]$index, 0)
  expect_length(nn_spherical(points, query, k = 99, search_radius = 180,
                             rule = "all", parallel = FALSE)[[1]]$index, 0)
  expect_identical(nn_spherical(points, query, k = 99,
                                parallel = FALSE)[[1]]$index, 1:5)
  legacy <- points
  names(legacy) <- c("zen", "azi")
  expect_identical(nn_spherical(legacy, legacy[1, ], k = NULL,
                                search_radius = 20, parallel = FALSE,
                                angle_columns = c("zen", "azi"))[[1]], all)
})

test_that("poles, missing directions and empty inputs retain row correspondence", {
  points <- data.frame(zenith = c(NA, 0, pi / 4, pi, pi / 2),
                       azimuth = c(0, NA, NA, NA, 0))
  neighbors <- nn_spherical(points, points, k = NULL, search_radius = 180,
                             parallel = FALSE)
  expect_identical(neighbors[[1]], list(index = integer(), distance = numeric()))
  expect_identical(neighbors[[3]], neighbors[[1]])
  expect_identical(neighbors[[2]]$index, c(2L, 5L, 4L))
  expect_equal(neighbors[[2]]$distance, c(0, pi / 2, pi))
  expect_identical(neighbors[[4]]$index, c(4L, 5L, 2L))
  expect_identical(nn_spherical(points, points[FALSE, ], parallel = FALSE), list())
  empty <- nn_spherical(points[FALSE, ], points, parallel = FALSE)
  expect_length(empty, nrow(points))
  expect_true(all(lengths(lapply(empty, function(x) x$index)) == 0))
})

test_that("invalid units and unbounded radius requests are rejected", {
  points <- data.frame(zenith = 0, azimuth = 0)
  expect_error(nn_spherical(points, points, k = NULL), "search_radius")
  expect_error(nn_spherical(points, points, search_radius = 181), "180 degrees")
  expect_error(nn_spherical(points, points, k = 0), "positive")
  bad <- data.frame(zenith = 45, azimuth = 0)
  expect_error(nn_spherical(bad, points), "radians")
  bad <- data.frame(zenith = 0, azimuth = 2 * pi)
  expect_error(nn_spherical(points, bad), "radians")
  expect_error(nn_spherical(points, points,
                            angle_columns = c("off_axis", "internal_azimuth")),
               "off_axis")
})
