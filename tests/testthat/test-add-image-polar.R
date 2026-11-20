test_that("image polar predictors use continuous image axes and dimensions", {
  sampling_coord <- data.frame(
    x = c(8, 4, 0, 4, 4, 8, 8, 10),
    y_down = c(3, 6, 3, 0, 3, 6, 3, 3),
    id = letters[1:8],
    zen = seq(0.1, 0.8, by = 0.1),
    azi = seq(0.1, 0.8, by = 0.1)
  )
  rownames(sampling_coord) <- rev(letters[1:8])

  observed <- add_image_polar(sampling_coord, width = 8, height = 6)

  expect_identical(observed[names(sampling_coord)], sampling_coord)
  expect_equal(
    observed$relative_radius,
    c(0.8, 0.6, 0.8, 0.6, 0, 1, 0.8, 1.2)
  )
  expect_equal(
    observed$image_angle,
    c(0, pi / 2, pi, -pi / 2, 0, atan2(3, 4), 0, 0)
  )
})

test_that("declared optical centers share polar predictors with internal geometry", {
  points <- data.frame(x = c(5, 3, 1, 3, 3), y_down = c(2, 4, 2, 0, 2),
                       id = letters[1:5])
  polar <- add_image_polar(points, 8, 6, center_of_vision = c(3, 2))
  expect_equal(polar$relative_radius, c(rep(.4, 4), 0))
  expect_equal(polar$image_angle, c(0, pi/2, pi, -pi/2, 0))
  expect_identical(polar[names(points)], points)
  expect_equal(add_image_polar(points, 8, 6, c(4, 3)), add_image_polar(points, 8, 6))
  expect_error(add_image_polar(points, 8, 6, c(3, NA)), "center_of_vision")
  expect_equal(nrow(add_image_polar(points[FALSE, ], 8, 6, c(3, 2))), 0)
  p <- list(off_axis_coef = c(1.1, .1), sin_internal_azimuth_coef = c(0, 1),
            cos_internal_azimuth_coef = c(-1, 0), atan_scale = .6)
  g <- set_internal_geometry(c(8, 6), list(model = "anisotropic_projection",
    center_of_vision = c(3, 2), parameters = p))
  coord <- coord_from_indices(c(1, 3, 5, 8), c(1, 2, 4, 6))
  prepared <- add_image_polar(coord, 8, 6, c(3, 2))
  predicted <- add_directions(indices_from_coord(coord$x, coord$y_down), g)
  t <- atan(prepared$relative_radius / .6)
  expect_equal(predicted$off_axis, 1.1 * t + .1 * t^2)
  expect_equal(predicted$internal_azimuth,
                atan2(cos(prepared$image_angle), -sin(prepared$image_angle)) %% (2*pi))
})

test_that("cell-center conversion and proportional scaling are explicit", {
  sampling_coord <- coord_from_indices(col = c(1, 8), row = c(1, 6))
  observed <- add_image_polar(sampling_coord, 8, 6)

  expect_equal(observed$relative_radius, rep(sqrt(18.5) / 5, 2))
  expect_equal(observed$image_angle, atan2(c(-2.5, 2.5), c(-3.5, 3.5)))
  expect_identical(observed[c("x", "y_down")], sampling_coord)

  scaled_coord <- sampling_coord * 2
  scaled <- add_image_polar(scaled_coord, 16, 12)
  expect_equal(
    scaled[c("relative_radius", "image_angle")],
    observed[c("relative_radius", "image_angle")]
  )
})

test_that("image preprocessing validates its coord contract", {
  sampling_coord <- data.frame(x = 1.5, y_down = 2.5)
  expect_error(add_image_polar(data.frame(col = 1, row = 2), 8, 6), "x, y_down")
  expect_error(add_image_polar(sampling_coord, 8.5, 6), "integer-like")
  expect_error(add_image_polar(sampling_coord, 8, 0), "positive")
  expect_error(
    add_image_polar(data.frame(x = Inf, y_down = 2), 8, 6),
    "Inf"
  )
  expect_error(
    add_image_polar(data.frame(x = NA_real_, y_down = 2), 8, 6),
    "NA"
  )
  prepared <- add_image_polar(sampling_coord, 8, 6)
  expect_error(add_image_polar(prepared, 8, 6), "reserved columns")
})
