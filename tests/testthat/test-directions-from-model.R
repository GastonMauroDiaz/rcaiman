# These tests protect the declared model and shared preprocessing, not the
# scientific validity of a camera calibration. The paired-image experiment
# still requires human review (design recorded on 2026-09-04).

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


test_that("model directions preserve rows and the declared sine convention", {
  sampling_coord <- add_image_polar(
    data.frame(
      x = c(8, 4, 0, 4, 8),
      y_down = c(3, 6, 3, 0, 3),
      dn = c(0.2, 0.4, 0.6, 0.8, 0.2),
      id = c("right", "bottom", "left", "top", "right again")
    ),
    8, 6
  )
  # These coefficients encode the package's North-up viewing convention.
  parameters <- list(
    zen_coef = 1,
    azi_sin_coef = c(0, -1),
    azi_cos_coef = c(-1, 0),
    atan_scale = 0.3
  )
  observed <- add_directions_from_model(
    sampling_coord,
    parameters
  )

  expect_identical(observed[names(sampling_coord)], sampling_coord)
  expect_equal(observed$zen, atan(c(0.8, 0.6, 0.8, 0.6, 0.8) / 0.3))
  expect_equal(observed$azi, c(3 * pi / 2, pi, pi / 2, 0, 3 * pi / 2))
  expect_true(all(observed$azi >= 0 & observed$azi < 2 * pi))

  target_coord <- reproject_points(observed, equidistant_coef, 100)
  expect_true(all(is.finite(target_coord$x)))
  expect_true(all(is.finite(target_coord$y_down)))
  expect_identical(target_coord$dn, observed$dn)
  expect_identical(target_coord$id, observed$id)
  expect_equal(target_coord$zen, observed$zen)
  expect_equal(target_coord$azi, observed$azi)
  expect_false(any(
    c("relative_radius", "image_angle") %in% names(target_coord)
  ))
})


test_that("shared preprocessing reproduces fitted model predictions", {
  # Synthetic observations test fit/evaluation consistency, not calibration.
  calibration_coord <- data.frame(
    x = c(1, 2, 3, 5, 7, 6, 8, 2, 1, 5, 7, 4),
    y_down = c(1, 5, 2, 4, 1, 5, 3, 2, 4, 1, 6, 5),
    zen = c(1.1, 0.8, 0.2, 0.4, 0.9, 0.7, 1.2, 0.5, 1, 0.6, 1.3, 0.3),
    azi = c(0.1, 1.2, 2.3, 3.4, 4.5, 5.6, 0.8, 1.9, 3, 4.1, 5.2, 6.1)
  )
  calibration_coord <- add_image_polar(calibration_coord, 8, 6)
  atan_scale <- 0.45
  calibration_coord$atan_radius <-
    atan(calibration_coord$relative_radius / atan_scale)
  zen_fit <- stats::lm(
    zen ~ poly(atan_radius, 3, raw = TRUE) - 1,
    data = calibration_coord
  )
  sin_fit <- stats::lm(
    sin(azi) ~ sin(image_angle) + cos(image_angle) - 1,
    data = calibration_coord
  )
  cos_fit <- stats::lm(
    cos(azi) ~ sin(image_angle) + cos(image_angle) - 1,
    data = calibration_coord
  )

  prediction_coord <- add_image_polar(
    data.frame(x = c(1.5, 4.5, 6.2), y_down = c(2.8, 1.7, 4.1)),
    8, 6
  )
  parameters <- list(
    zen_coef = stats::coef(zen_fit),
    azi_sin_coef = stats::coef(sin_fit),
    azi_cos_coef = stats::coef(cos_fit),
    atan_scale = atan_scale
  )
  observed <- add_directions_from_model(prediction_coord, parameters)
  prediction_coord$atan_radius <-
    atan(prediction_coord$relative_radius / atan_scale)

  expect_equal(
    observed$zen,
    unname(stats::predict(zen_fit, prediction_coord))
  )
  expect_equal(
    observed$azi,
    unname(atan2(
      stats::predict(sin_fit, prediction_coord),
      stats::predict(cos_fit, prediction_coord)
    ) %% (2 * pi))
  )
})


test_that("image model functions accept structured empty inputs", {
  sampling_coord <- data.frame(
    x = numeric(), y_down = numeric(), id = character()
  )
  prepared <- add_image_polar(sampling_coord, 8, 6)
  parameters <- list(
    zen_coef = 1,
    azi_sin_coef = c(0, -1),
    azi_cos_coef = c(-1, 0),
    atan_scale = 0.6
  )
  observed <- add_directions_from_model(prepared, parameters)

  expect_identical(observed[names(sampling_coord)], sampling_coord)
  expect_identical(observed$relative_radius, numeric())
  expect_identical(observed$image_angle, numeric())
  expect_identical(observed$zen, numeric())
  expect_identical(observed$azi, numeric())
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


test_that("model evaluation validates coefficients and preserves raw zenith", {
  prepared <- add_image_polar(data.frame(x = 1.5, y_down = 2.5), 8, 6)
  parameters <- list(
    zen_coef = 1,
    azi_sin_coef = c(0, -1),
    azi_cos_coef = c(-1, 0),
    atan_scale = 0.6
  )
  args <- list(sampling_coord = prepared, parameters = parameters)
  parameters$atan_scale <- 0
  expect_error(
    .assert_anisotropic_projection_parameters(parameters),
    "positive"
  )
  parameters$atan_scale <- Inf
  expect_error(
    .assert_anisotropic_projection_parameters(parameters),
    "Inf"
  )
  parameters$atan_scale <- 0.6
  parameters$zen_coef <- numeric()
  expect_error(.assert_anisotropic_projection_parameters(parameters), "empty")
  parameters$zen_coef <- c(1, 0, Inf)
  expect_error(.assert_anisotropic_projection_parameters(parameters), "Inf")
  parameters$zen_coef <- c(-0.1, 0, 0, 0, 0)
  args$parameters <- parameters
  observed <- do.call(add_directions_from_model, args)
  expect_lt(observed$zen, 0)
  args$sampling_coord <- observed
  expect_error(do.call(add_directions_from_model, args), "reserved columns")
  args$sampling_coord <- prepared[c("x", "y_down")]
  expect_error(do.call(add_directions_from_model, args), "relative_radius")
})


test_that("anisotropic projection contract is named and has no zenith intercept", {
  parameters <- list(
    zen_coef = c(1, 2, 3, 4, 5),
    azi_sin_coef = c(0, -1),
    azi_cos_coef = c(-1, 0),
    atan_scale = 0.6
  )
  expect_true(.assert_anisotropic_projection_parameters(parameters))
  expect_true(
    .assert_anisotropic_projection_parameters(parameters[rev(names(parameters))])
  )

  at_center <- add_image_polar(
    data.frame(x = 4, y_down = 3), width = 8, height = 6
  )
  observed <- add_directions_from_model(at_center, parameters)
  expect_equal(observed$zen, 0)

  expect_error(
    .assert_anisotropic_projection_parameters(unname(parameters)),
    "named list"
  )
  parameters$unexpected <- 1
  expect_error(.assert_anisotropic_projection_parameters(parameters), "exactly")
})
