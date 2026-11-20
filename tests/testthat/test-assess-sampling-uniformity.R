sampling_area_fixture <- function() {
  seg <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 3))
  seg[] <- c(10, 10, 90, 90, 0, NA)
  w <- terra::setValues(terra::rast(seg), c(1, 1, 3, 3, NA, NA))
  list(seg = seg, w = w, indices = data.frame(col = c(1, 3), row = c(1, 1)))
}

test_that("angular areas define expected counts and coverage, without metadata", {
  f <- sampling_area_fixture()
  # Segment areas are 2 and 6, despite arbitrary IDs and raster name.
  balanced <- f$indices[c(1, 2, 2, 2), ]
  expect_equal(assess_sampling_uniformity(balanced, f$seg, f$w),
               list(kl = 0, coverage = 1))
  expect_equal(assess_sampling_uniformity(f$indices, f$seg, f$w),
               list(kl = .5 * log(2) + .5 * log(2/3), coverage = 1))
  expect_equal(assess_sampling_uniformity(f$indices[1, ], f$seg, f$w),
               list(kl = log(4), coverage = .25))
  expect_equal(assess_sampling_uniformity(balanced, f$seg, f$w * 1e200),
               assess_sampling_uniformity(balanced, f$seg, f$w))
  expect_equal(assess_sampling_uniformity(balanced[FALSE, ], f$seg, f$w),
               list(kl = NA_real_, coverage = 0))
  balanced$segment <- 123 # Additional columns do not collide with raster names.
  expect_equal(assess_sampling_uniformity(balanced, f$seg, f$w)$kl, 0)
})

test_that("an explicit domain restriction cuts segments but gaps do not change it", {
  f <- sampling_area_fixture()
  m <- terra::setValues(terra::rast(f$seg), c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(assess_sampling_uniformity(f$indices[1, ], f$seg, f$w, m),
               list(kl = log(5/2), coverage = 2/5))
  f$seg[4] <- 0
  expect_equal(assess_sampling_uniformity(f$indices, f$seg, f$w),
               assess_sampling_uniformity(f$indices, f$seg, f$w, m))
  m[3] <- FALSE
  expect_equal(assess_sampling_uniformity(f$indices[1, ], f$seg, f$w, m),
               list(kl = 0, coverage = 1))
  expect_error(assess_sampling_uniformity(f$indices, f$seg, f$w, m), "reference domain")
})

test_that("invalid indices, weights, masks and grids are rejected", {
  f <- sampling_area_fixture()
  for (ind in list(data.frame(col = 2, row = 2), data.frame(col = 4, row = 1),
                   data.frame(col = 3, row = 2))) {
    expect_error(assess_sampling_uniformity(ind, f$seg, f$w), "positive cell IDs")
  }
  for (value in c(-1, Inf)) {
    bad <- terra::deepcopy(f$w); bad[1] <- value
    expect_error(assess_sampling_uniformity(f$indices, f$seg, bad), "w")
  }
  for (value in c(0, NA)) {
    bad <- terra::deepcopy(f$w); bad[1] <- value
    expect_error(assess_sampling_uniformity(f$indices, f$seg, bad), "reference domain")
  }
  expect_error(assess_sampling_uniformity(f$indices, f$seg, f$w * 0), "positive weights")
  m <- terra::setValues(terra::rast(f$seg), TRUE); m[1] <- NA
  expect_error(assess_sampling_uniformity(f$indices, f$seg, f$w, m), "NA-free")
  bad <- terra::deepcopy(f$w); terra::ext(bad) <- terra::ext(0, 30, 0, 20)
  expect_error(assess_sampling_uniformity(f$indices, f$seg, bad), "same geometry")
})

test_that("Beckers and HEALPix use actual areas for tilted and rectangular FOVs", {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(40, 30), normalization_radius = 30,
    normalization_angle = 90)
  for (dims in list(c(80, 60), c(60, 40))) {
    g <- set_internal_geometry(dims, geometry)
    w <- area_weight_image(g, geometry)
    for (elevation in c(90, 32.5)) {
      d <- direction_image(g, elevation = elevation)
      for (seg in list(beckers_segmentation(d, 20), healpix_segmentation(d, 2))) {
        labels <- terra::values(seg, mat = FALSE)
        weights <- terra::values(w, mat = FALSE)
        selected <- which(labels > 0 & !is.na(weights))
        areas <- tapply(weights[selected], labels[selected], sum)
        q <- areas / sum(areas)
        cells <- selected[!duplicated(labels[selected])]
        rc <- terra::rowColFromCell(seg, cells)
        indices <- data.frame(col = rc[, 2], row = rc[, 1])
        result <- assess_sampling_uniformity(indices, seg, w)
        expect_equal(result$coverage, 1)
        expect_equal(result$kl, mean(log((1/length(q)) / q)))
        # Raster names and constructor metadata have no role in the metrics.
        names(seg) <- "arbitrary"
        attr(seg, "n_cells") <- NULL
        attr(seg, "cell_solid_angle") <- NULL
        expect_equal(assess_sampling_uniformity(indices, seg, w), result)
      }
    }
  }
})

test_that("uniform solid-angle sampling matches clipped HEALPix reference", {
  # Equal steps in cos(zenith) and azimuth give equal-solid-angle samples.
  d <- set_rcaiman_geometry(terra::rast(nrows = 60, ncols = 120, nlyrs = 2))
  d[] <- cbind(rep(acos(((1:60)-.5)/60)*180/pi, each = 120),
                rep(((1:120)-.5)*360/120, 60))
  names(d) <- c("zenith", "azimuth")
  seg <- healpix_segmentation(d, 1)
  w <- terra::setValues(terra::rast(seg), 1)
  indices <- expand.grid(col = 1:120, row = 1:60)
  expect_equal(assess_sampling_uniformity(indices, seg, w), list(kl = 0, coverage = 1))
  q <- prop.table(table(terra::values(seg, mat = FALSE)))
  expect_equal(as.numeric(q[1:4]), rep(1/6, 4), tolerance = .001)
  expect_equal(as.numeric(q[5:8]), rep(1/12, 4), tolerance = .001)
  expect_gt(sum(q * log(q * length(q))), .05)
})


test_that("equal represented areas recover equal cell probabilities", {
  f <- sampling_area_fixture()
  f$w[] <- 1
  expect_equal(assess_sampling_uniformity(f$indices, f$seg, f$w),
               list(kl = 0, coverage = 1))
  expect_equal(assess_sampling_uniformity(f$indices[1, ], f$seg, f$w),
               list(kl = log(2), coverage = .5))
})
