sun_estimation_fixture <- function(size = 120, n_cells = 400, elevation = 90) {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = rep(size/2, 2), normalization_radius = size/2 - 2,
    normalization_angle = 90)
  g <- set_internal_geometry(rep(size, 2), geometry)
  d <- direction_image(g, elevation = elevation)
  w <- area_weight_image(g, geometry)
  seg <- beckers_segmentation(d, n_cells)
  domain <- seg > 0
  list(d = d, w = w, seg = seg, bin = domain,
       r = terra::setValues(terra::rast(w), 1))
}

test_that("veiled mode omits vegetation and uses the entire cell's center", {
  f <- sun_estimation_fixture()
  labels <- terra::values(f$seg, mat = FALSE)
  first <- which(labels == 2)
  second <- which(labels == 3)
  f$r[] <- 0
  f$r[first] <- 1
  f$r[tail(first, 1)] <- 100
  f$r[second] <- 80
  f$bin[] <- FALSE
  f$bin[c(tail(first, 1), second)] <- TRUE
  before <- terra::values(f$r)
  out <- estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w,
                            mode = "veiled", return = "all")
  cell <- out$cells[out$cells$id == 2, ]
  expect_equal(cell$brightness, 100)
  expect_equal(cell$n_samples, 1)
  expect_identical(out$method, "brightest_cell")
  expect_equal(unname(as.numeric(out$sun_angles[1, ])), c(cell$zenith, cell$azimuth))
  expected <- compute_features(f$r, f$seg, f$d, f$w)
  i <- match(2, expected$id)
  expect_equal(cell$area, expected$area[i])
  expect_equal(cell$zenith, expected$spherical_centroid_zenith[i])
  expect_equal(cell$azimuth, expected$spherical_centroid_azimuth[i])
  expect_equal(terra::values(f$r), before)
  expect_equal(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, mode = "veiled"),
               out$sun_angles)
})

test_that("aggregation and cell geometry are independent of within-cell gap location", {
  f <- sun_estimation_fixture()
  labels <- terra::values(f$seg, mat = FALSE)
  # Constant measurements per cell isolate the geometry from radiometric changes.
  f$r[] <- ifelse(labels > 0, labels, NA)
  full <- estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, return = "all")
  ids <- which(labels > 0)
  for (last in c(FALSE, TRUE)) {
    cells <- ids[!duplicated(labels[ids], fromLast = last)]
    f$bin[] <- FALSE; f$bin[cells] <- TRUE
    sparse <- estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, return = "all")
    fields <- setdiff(names(full$cells), "n_samples")
    expect_equal(sparse$cells[fields], full$cells[fields])
    expect_equal(sparse$groups, full$groups)
    expect_equal(sparse$sun_angles, full$sun_angles)
    expect_equal(sparse$centroid, full$centroid)
  }
})

test_that("no samples, invalid inputs, and uniform brightness have explicit outcomes", {
  f <- sun_estimation_fixture()
  out <- estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, return = "all")
  expect_identical(out$method, "spherical_centroid")
  expect_false("circle" %in% names(out))
  expect_true(all(out$cells$bright))
  expect_equal(out$sun_angles, out$centroid)
  f$bin[] <- FALSE
  expect_error(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w), "No sky measurements")
  f$bin <- f$seg > 0; f$r[] <- NA_real_
  expect_error(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w), "No sky measurements")
  f$r[] <- Inf
  expect_error(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w), "sky measurements")
  f$r[] <- 1
  expect_error(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, 181), "circumsolar_angle")
  expect_error(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, mode = "bad"), "mode")
  f$bin[1] <- NA
  expect_error(estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w), "NA-free")
})

test_that("area replaces cell count when ranking equally bright groups", {
  d <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 7, nlyrs = 2))
  d[] <- cbind(rep(40, 7), seq(10, 70, 10)); names(d) <- c("zenith", "azimuth")
  seg <- terra::setValues(terra::rast(d[[1]]), c(1, 2, 0, 3, 0, 4, 5))
  r <- terra::setValues(terra::rast(seg), c(100, 100, 0, 100, 0, 1, 1))
  w <- terra::setValues(terra::rast(seg), c(1, 1, 1, 5, 1, 1, 1))
  out <- estimate_sun_angles(r, d, seg > 0, seg, w, circumsolar_angle = 1, return = "all")
  expect_equal(out$groups$area, c(2, 5))
  expect_equal(out$groups$seed, c(FALSE, TRUE))
  expect_equal(out$cells$id[out$cells$selected], 3)
  expect_identical(out$method, "spherical_centroid")
  expect_equal(out$sun_angles$azimuth, 40*pi/180)
})

test_that("domain restrictions change cell areas independently of canopy gaps", {
  f <- sun_estimation_fixture(elevation = 32.5)
  f$r[] <- seq_len(terra::ncell(f$r))
  base <- estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, return = "all")
  m <- !is.na(f$d$zenith) & f$d$zenith <= 70
  out <- estimate_sun_angles(f$r, f$d, f$bin, f$seg, f$w, m = m, return = "all")
  expect_lt(sum(out$cells$area), sum(base$cells$area))
  expect_true(all(out$cells$zenith < 70*pi/180))
  expect_false(any(c("boundary", "contour", "hull", "resolution") %in% names(out$cells)))
  # Domain area uses all represented pixels, not just the sky measurements.
  labels <- terra::values(f$seg, mat = FALSE)
  weights <- terra::values(f$w, mat = FALSE)
  domain <- terra::values(m, mat = FALSE) & !is.na(labels) & labels > 0 &
    !is.na(weights) & weights > 0
  expect_equal(sum(out$cells$area), sum(weights[domain]))
})

test_that("CIE twilight skies return the visible centroid with either partition", {
  f <- sun_estimation_fixture(size = 320, n_cells = 2400)
  truth <- data.frame(zenith = 92*pi/180, azimuth = 0)
  f$r <- cie_image(f$d, truth, as.numeric(cie_table[10, 1:5]))
  cb <- chessboard(f$r, 10)
  set.seed(123)
  f$bin <- as.logical(max(cb == sample(unique(cb[]), 300)))
  for (seg in list(f$seg, healpix_segmentation(f$d, 16))) {
    out <- estimate_sun_angles(f$r, f$d, f$bin, seg, f$w, return = "all")
    expect_identical(out$method, "spherical_centroid")
    expect_equal(out$sun_angles, out$centroid)
    expect_lt(out$sun_angles$zenith, pi/2)
    expect_identical(names(out), c("sun_angles", "method", "centroid", "cells", "groups"))
    expect_identical(names(out$cells), c("id", "area", "zenith", "azimuth",
      "n_samples", "brightness", "group", "bright", "selected"))
    expect_true(all(out$cells$selected <= out$cells$bright))
    # Independently aggregate every represented pixel of the selected cells.
    # This includes pixels hidden by the canopy, as the contract requires.
    labels <- terra::values(seg, mat = FALSE)
    selected <- labels %in% out$cells$id[out$cells$selected]
    z <- terra::values(f$d$zenith, mat = FALSE)[selected] * pi/180
    a <- terra::values(f$d$azimuth, mat = FALSE)[selected] * pi/180
    weights <- terra::values(f$w, mat = FALSE)[selected]
    x <- sum(weights * sin(z) * cos(a))
    y <- sum(weights * sin(z) * sin(a))
    up <- sum(weights * cos(z))
    expected <- data.frame(zenith = atan2(sqrt(x*x + y*y), up),
                           azimuth = atan2(y, x) %% (2*pi))
    expect_equal(out$sun_angles, expected, tolerance = 1e-12)
    expect_true(any(!terra::values(f$bin, mat = FALSE)[selected]))
    expect_equal(estimate_sun_angles(f$r, f$d, f$bin, seg, f$w), out$sun_angles)
    expect_identical(names(out$sun_angles), c("zenith", "azimuth"))
    expect_equal(nrow(out$sun_angles), 1)
  }
})

test_that("zero resultant directions cannot invent a cell center", {
  d <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = 2, nlyrs = 2))
  d[] <- cbind(c(90, 90), c(0, 180)); names(d) <- c("zenith", "azimuth")
  r <- terra::setValues(terra::rast(d[[1]]), 1)
  expect_error(estimate_sun_angles(r, d, r > 0, r, r), "undefined spherical centroid")
})
