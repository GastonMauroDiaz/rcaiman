cone_samples <- function() {
  x <- sky_grid_centers(15)
  x$Blue <- with(x, 100 + 20 * zenith - 3 * zenith^2 +
                  4 * sin(azimuth) - 7 * cos(azimuth)) +
    sin(seq_len(nrow(x))) / 10
  x
}

cone_directions <- function() {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
                   center_of_vision = c(20, 15), normalization_radius = 15,
                   normalization_angle = 90)
  g <- set_internal_geometry(c(40, 30), geometry)
  direction_image(g, sensor_rotation = 12, azimuth_rotation = 31, elevation = 40)
}

test_that("fitting preserves both formulas with external angles in radians", {
  x <- cone_samples()
  for (method in c("zenith_only", "zenith_n_azimuth")) {
    fit <- fit_cone_shaped_model(x, "Blue", method, return = "all")
    old <- if (method == "zenith_only") {
      lm(Blue ~ poly(zenith, 2, raw = TRUE), data = x)
    } else {
      lm(Blue ~ poly(zenith, 2, raw = TRUE) + sin(azimuth) + cos(azimuth), data = x)
    }
    expect_named(fit, c("parameters", "model"))
    expect_s3_class(fit$model, "lm")
    expect_identical(class(fit$parameters), "list")
    expect_identical(fit$parameters$method, method)
    expect_equal(fit$parameters$coefficients, coef(fit$model))
    expect_equal(unname(coef(fit$model)), unname(coef(old)))
    expect_equal(fitted(fit$model), fitted(old))
    expect_equal(residuals(fit$model), residuals(old))
    expect_equal(predict(fit$model, newdata = x), predict(old, newdata = x))
    expect_equal(fit_cone_shaped_model(x, 3L, method), fit$parameters)
  }
  # No indices, or even an azimuth column, are needed for zenith-only fitting.
  expect_equal(fit_cone_shaped_model(x[c("zenith", "Blue")], "Blue", "zenith_only"),
               fit_cone_shaped_model(x, "Blue", "zenith_only"))
})

test_that("insufficient and nonidentifiable fits warn and return NULL", {
  x <- cone_samples()
  expect_warning(result <- fit_cone_shaped_model(x[1:19, ], "Blue"), "20")
  expect_null(result)
  x$zenith <- 1
  expect_warning(result <- fit_cone_shaped_model(x, "Blue"), "full-rank")
  expect_null(result)
  expect_warning(result <- fit_cone_shaped_model(x, "Blue", "zenith_only"), "full-rank")
  expect_null(result)
  x <- cone_samples()
  x$azimuth <- 0
  expect_warning(result <- fit_cone_shaped_model(x, "Blue"), "full-rank")
  expect_null(result)
})

test_that("fitting requires finite observations and appropriately named angles", {
  x <- cone_samples()
  expect_error(fit_cone_shaped_model(x, 2.5), "integer")
  expect_error(fit_cone_shaped_model(x, "missing"), "required columns")
  expect_error(fit_cone_shaped_model(x, "Blue", return = "fun"))
  expect_error(fit_cone_shaped_model(x, "Blue", method = "other"))
  missing <- x
  names(missing)[1] <- "zen"
  expect_error(fit_cone_shaped_model(missing, "Blue"), "zenith")
  x$Blue[1] <- NA_real_
  expect_error(fit_cone_shaped_model(x, "Blue"), "NA")
  x <- cone_samples()
  x$zenith[1] <- 90
  expect_error(fit_cone_shaped_model(x, "Blue"), "radians")
  x$zenith[1] <- 0
  x$azimuth[1] <- NA_real_
  expect_error(fit_cone_shaped_model(x, "Blue"), "NA")
  expect_type(fit_cone_shaped_model(x, "Blue", "zenith_only"), "list")
  x$azimuth[1] <- 2 * pi
  expect_error(fit_cone_shaped_model(x, "Blue"), "radians")
  x$azimuth[1] <- Inf
  expect_error(fit_cone_shaped_model(x, "Blue"), "Inf")
})

test_that("raster predictions agree with lm on an inclined acquisition", {
  d <- cone_directions()
  query <- as.data.frame(terra::values(d)) * pi / 180
  m <- !is.na(d$zenith) & d$zenith < 80
  for (method in c("zenith_only", "zenith_n_azimuth")) {
    fit <- fit_cone_shaped_model(cone_samples(), "Blue", method, return = "all")
    expected <- as.numeric(predict(fit$model, newdata = query))
    out <- cone_shaped_image(d, fit$parameters)
    expect_equal(terra::values(out, mat = FALSE), expected)
    expect_true(terra::compareGeom(out, d[[1]]))
    expect_identical(names(out), "prediction")
    expected[!terra::values(m, mat = FALSE)] <- NA_real_
    masked <- cone_shaped_image(d, fit$parameters, m)
    expect_equal(terra::values(masked, mat = FALSE), expected)
    reordered <- fit$parameters
    reordered$coefficients <- rev(reordered$coefficients)
    expect_equal(terra::values(cone_shaped_image(d, reordered)),
                 terra::values(out))
  }
})

test_that("the sampling and fitting workflow reconstructs a known sky", {
  d <- cone_directions()
  r <- 100 + 20 * d$zenith * pi / 180 - 3 * (d$zenith * pi / 180)^2 +
    4 * sin(d$azimuth * pi / 180) - 7 * cos(d$azimuth * pi / 180)
  names(r) <- "Blue"
  m <- !is.na(d$zenith) & d$zenith < 90
  indices <- expand.grid(col = seq(2, 38, 3), row = seq(2, 29, 3))
  indices <- rem_unselected(indices, m)
  samples <- add_directions(extract_dn(r, indices, use_window = FALSE), d = d)
  parameters <- fit_cone_shaped_model(samples, "Blue")
  expect_equal(unname(parameters$coefficients), c(100, 20, -3, 4, -7))
  reconstructed <- cone_shaped_image(d, parameters)
  expect_equal(terra::values(reconstructed), terra::values(r), ignore_attr = TRUE)
})

test_that("missing directions and masks follow the selected model", {
  d <- cone_directions()
  z <- rep(c(0, 45, 90, NA_real_), length.out = terra::ncell(d))
  a <- rep(c(NA_real_, 180, 90, 0), length.out = terra::ncell(d))
  terra::values(d) <- cbind(z, a)
  for (method in c("zenith_only", "zenith_n_azimuth")) {
    parameters <- fit_cone_shaped_model(cone_samples(), "Blue", method)
    values <- terra::values(cone_shaped_image(d, parameters), mat = FALSE)
    expect_identical(is.na(values), if (method == "zenith_only") is.na(z) else
                       is.na(z) | is.na(a))
    m <- !is.na(d$zenith) & d$zenith < 0
    expect_true(all(is.na(terra::values(cone_shaped_image(d, parameters, m)))))
  }
  parameters <- fit_cone_shaped_model(cone_samples(), "Blue", "zenith_only")
  parameters$coefficients[] <- c(-5, 0, 0)
  out <- terra::values(cone_shaped_image(d, parameters), mat = FALSE)
  expect_true(all(out[!is.na(z)] == -5))
  terra::values(d) <- NA_real_
  expect_true(all(is.na(terra::values(cone_shaped_image(d, parameters)))))
})

test_that("prediction rejects malformed parameters and incompatible masks", {
  d <- cone_directions()
  parameters <- fit_cone_shaped_model(cone_samples(), "Blue")
  expect_error(cone_shaped_image(d, list()), "named list")
  bad <- parameters
  bad$coefficients[1] <- NA_real_
  expect_error(cone_shaped_image(d, bad), "NA")
  bad$coefficients <- unname(parameters$coefficients)
  expect_error(cone_shaped_image(d, bad), "names")
  bad$coefficients <- parameters$coefficients[-1]
  expect_error(cone_shaped_image(d, bad), "length")
  expect_error(cone_shaped_image(d[[1]], parameters), "two layers")
  m <- !is.na(d$zenith)
  terra::values(m)[1] <- NA
  expect_error(cone_shaped_image(d, parameters, m))
  m <- terra::rast(nrows = 2, ncols = 2)
  terra::values(m) <- TRUE
  expect_error(cone_shaped_image(d, parameters, m))
})
