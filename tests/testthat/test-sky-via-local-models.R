local_sky_fixture <- function() {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(16, 12), normalization_radius = 12,
    normalization_angle = 90)
  g <- set_internal_geometry(c(32, 24), geometry)
  d <- direction_image(g, azimuth_rotation = 37, elevation = 60)
  samples <- as.data.frame(.cartesian_to_directions(.fibonacci(7)))
  names(samples) <- c("zenith", "azimuth")
  samples$dn <- 1000 + 200 * samples$zenith^2 +
    80 * sin(samples$azimuth) + 50 * cos(2 * samples$azimuth)
  list(d = d, samples = samples)
}

test_that("local reconstruction uses external directions on a tilted rectangular grid", {
  f <- local_sky_fixture()
  samples <- f$samples
  samples$dn <- 7
  sky <- sky_via_local_models(f$d, samples, "dn", spacing_deg = 20, spacing_px = 1)
  selected <- with(as.data.frame(f$d, na.rm = FALSE),
    !is.na(zenith) & !is.na(azimuth) & zenith <= 90)
  v <- terra::values(sky, mat = FALSE)
  expect_true(terra::compareGeom(sky, f$d[[1]]))
  expect_identical(names(sky), "sky")
  expect_equal(v[selected], rep(7, sum(selected)), tolerance = 1e-9)
  expect_true(all(is.na(v[!selected])))
  samples$col <- NA_integer_
  samples$row <- -1L
  same <- sky_via_local_models(f$d, samples, 3, spacing_deg = 20, spacing_px = 1)
  expect_equal(terra::values(same), terra::values(sky))
})

test_that("conicity one interpolates predictions of the global fit", {
  f <- local_sky_fixture()
  fit <- fit_cone_shaped_model(f$samples, "dn", return = "all")
  queries <- as.data.frame(.cartesian_to_directions(.fibonacci(20)))
  names(queries) <- c("zenith", "azimuth")
  queries$prediction <- as.numeric(predict(fit$model, newdata = queries))
  indices <- expand.grid(col = 1:32, row = 1:24)
  m <- !is.na(f$d$zenith) & !is.na(f$d$azimuth) & f$d$zenith <= 90
  indices <- rem_unselected(indices, m)
  directions <- add_directions(indices, d = f$d)
  expected <- interpolate_spherical(queries, directions, "prediction",
    k = 3, p = 2, search_radius = 40)
  sky <- sky_via_local_models(f$d, f$samples, "dn", conicity = 1,
    spacing_deg = 20, spacing_px = 1)
  actual <- extract_dn(sky, indices, use_window = FALSE)$sky
  expect_equal(actual, expected$prediction, tolerance = 1e-8)
  local <- sky_via_local_models(f$d, f$samples, "dn", conicity = 0,
    spacing_deg = 20, spacing_px = 1)
  expect_gt(max(abs(terra::values(local - sky)), na.rm = TRUE), 1)
})

test_that("mask limits reconstruction without filtering observations", {
  f <- local_sky_fixture()
  m <- !is.na(f$d$zenith) & f$d$zenith >= 20 & f$d$zenith <= 55
  masked <- sky_via_local_models(f$d, f$samples, "dn", conicity = .5,
    spacing_deg = 20, spacing_px = 1, m = m)
  full <- sky_via_local_models(f$d, f$samples, "dn", conicity = .5,
    spacing_deg = 20, spacing_px = 1)
  keep <- terra::values(m, mat = FALSE)
  expect_true(all(is.na(terra::values(masked, mat = FALSE)[!keep])))
  expect_equal(terra::values(masked, mat = FALSE)[keep],
    terra::values(full, mat = FALSE)[keep])
  empty <- m & FALSE
  expect_true(all(is.na(terra::values(
    sky_via_local_models(f$d, f$samples, "dn", m = empty)))))
})

test_that("missing azimuth at zenith is not assigned an arbitrary prediction", {
  f <- local_sky_fixture()
  f$d$zenith[1] <- 0
  f$d$azimuth[1] <- NA_real_
  sky <- sky_via_local_models(f$d, f$samples, "dn", spacing_deg = 30, spacing_px = 1)
  expect_true(is.na(terra::values(sky)[1]))
})

test_that("unsupported fits return NA and invalid samples are not silently removed", {
  f <- local_sky_fixture()
  expect_warning(sky <- sky_via_local_models(f$d, f$samples[1:19, ], "dn"),
    "At least 20")
  expect_true(all(is.na(terra::values(sky))))
  repeated <- f$samples[rep(1, 20), ]
  suppressWarnings(expect_warning(
    sky <- sky_via_local_models(f$d, repeated, "dn", spacing_deg = 180),
    "full-rank"))
  expect_true(all(is.na(terra::values(sky))))
  for (column in c("zenith", "azimuth", "dn")) {
    bad <- f$samples
    bad[[column]][1] <- NA_real_
    expect_error(sky_via_local_models(f$d, bad, "dn"), "NA")
  }
  bad <- f$samples
  bad$zenith[1] <- 45
  expect_error(sky_via_local_models(f$d, bad, "dn"), "radians")
  expect_error(sky_via_local_models(f$d, f$samples, "dn", conicity = 2), "conicity")
  expect_error(sky_via_local_models(f$d, f$samples, "dn", spacing_deg = 181), "spacing_deg")
  expect_error(sky_via_local_models(f$d, f$samples, "dn", spacing_px = 1.5), "integer")
})
