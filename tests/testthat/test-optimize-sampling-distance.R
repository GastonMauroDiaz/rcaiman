sampling_distance_fixture <- function() {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(40, 30), normalization_radius = 30,
    normalization_angle = 90)
  g <- set_internal_geometry(c(80, 60), geometry)
  d <- direction_image(g, elevation = 32.5)
  seg <- beckers_segmentation(d, 20)
  indices <- expand.grid(col = seq(1, 80, 3), row = seq(1, 60, 3))
  indices <- rem_unselected(indices, seg > 0)
  list(seg = seg, w = area_weight_image(g, geometry), indices = add_directions(indices, g = g, d = d))
}

test_that("optimization uses current directions and retains original metrics", {
  fixture <- sampling_distance_fixture()
  distances <- seq(1, 45, 8)
  result <- optimize_sampling_distance(fixture$indices, fixture$seg, fixture$w, distances, 3)
  expect_true(is.finite(result$optimum))
  expect_true(result$optimum >= min(distances) && result$optimum <= max(distances))
  raw <- lapply(distances, function(distance) {
    retained <- rem_nearby_indices(fixture$indices, distance, space = "spherical")
    assess_sampling_uniformity(retained, fixture$seg, fixture$w)
  })
  expect_equal(result$data$kl_raw, vapply(raw, `[[`, numeric(1), "kl"))
  expect_equal(result$data$coverage_raw, vapply(raw, `[[`, numeric(1), "coverage"))
  smooth <- function(x) (c(x[1], head(x, -1)) + x + c(tail(x, -1), tail(x, 1))) / 3
  expect_equal(result$data$kl, normalize_min_max(smooth(result$data$kl_raw)))
  expect_equal(result$data$coverage, normalize_min_max(smooth(result$data$coverage_raw)))
  internal <- optimize_sampling_distance(fixture$indices, fixture$seg, fixture$w, distances, 3,
    angle_columns = c("off_axis", "internal_azimuth"))
  expect_equal(internal, result)
})

test_that("unchanged sampling warns and preserves raw results", {
  fixture <- sampling_distance_fixture()
  expect_warning(result <- optimize_sampling_distance(fixture$indices, fixture$seg, fixture$w,
    c(.001, .002), window_length = 3), "kl and coverage did not vary")
  expect_identical(result$optimum, NA_real_)
  expect_true(all(is.finite(result$data$kl_raw)))
  expect_true(all(is.finite(result$data$coverage_raw)))
  expect_true(all(is.na(result$data[, c("kl", "coverage", "diff")])))
})

test_that("constant coverage does not discard a changing KL curve", {
  testthat::local_mocked_bindings(
    .sampling_area_reference = function(...) NULL,
    rem_nearby_indices = function(sampling_indices, min_dist, ...) data.frame(distance = min_dist),
    .assess_sampling_reference = function(sampling_indices, ...) {
      if (is.null(sampling_indices$distance)) return(list(kl = 0, coverage = 1))
      list(kl = sampling_indices$distance, coverage = .8)
    }
  )
  expect_warning(result <- optimize_sampling_distance(data.frame(col = 1, row = 1),
    NULL, NULL, 1:3, window_length = 3), "coverage did not vary")
  expect_identical(result$optimum, NA_real_)
  expect_equal(result$data$kl_raw, 1:3)
  expect_equal(result$data$coverage_raw, rep(.8, 3))
  expect_equal(result$data$kl, c(0, .5, 1))
  expect_true(all(is.na(result$data$coverage)))
  expect_true(all(is.na(result$data$diff)))
})

test_that("the first crossing is interpolated and identical curves select the first distance", {
  kl <- c(0, 1, .2, .9)
  coverage <- c(1, 0, 1, 0)
  testthat::local_mocked_bindings(
    .sampling_area_reference = function(...) NULL,
    rem_nearby_indices = function(sampling_indices, min_dist, ...) data.frame(distance = min_dist),
    .assess_sampling_reference = function(sampling_indices, ...) {
      i <- if (is.null(sampling_indices$distance)) 1 else sampling_indices$distance
      list(kl = kl[i], coverage = coverage[i])
    }
  )
  indices <- data.frame(col = 1, row = 1)
  result <- optimize_sampling_distance(indices, NULL, NULL, 1:4)
  expect_equal(result$optimum, 1.5)
  coverage <- kl
  result <- optimize_sampling_distance(indices, NULL, NULL, 1:4)
  expect_equal(result$optimum, 1)
})

test_that("invalid candidates and samples fail before thinning", {
  fixture <- sampling_distance_fixture()
  for (distances in list(numeric(0), 0, 91, c(2, 1), c(1, 1))) {
    expect_error(optimize_sampling_distance(fixture$indices, fixture$seg, fixture$w, distances), "distances")
  }
  expect_error(optimize_sampling_distance(fixture$indices[FALSE, ], fixture$seg, fixture$w, 1:3),
               "at least one row")
  expect_error(optimize_sampling_distance(fixture$indices, fixture$seg, fixture$w, 1:3, 2),
               "window_length")
  invalid <- rbind(fixture$indices[1, ], fixture$indices[1, ])
  invalid$col[2] <- 81
  expect_error(optimize_sampling_distance(invalid, fixture$seg, fixture$w, 1:3), "positive cell IDs")
})

test_that("distance search uses the same boundary priority as final thinning", {
  fixture <- sampling_distance_fixture()
  distance <- distance_to_black(fixture$seg > 0)
  names(distance) <- "distance_to_black"
  x <- extract_dn(distance, fixture$indices, use_window = FALSE)
  distances <- seq(1, 45, 8)
  for (priority in c("higher", "lower")) {
    result <- optimize_sampling_distance(x, fixture$seg, fixture$w, distances,
      col_id = "distance_to_black", prioritize = priority)
    raw <- lapply(distances, function(radius) {
      selected <- rem_nearby_indices(x, radius, col_id = "distance_to_black",
        prioritize = priority, space = "spherical")
      assess_sampling_uniformity(selected, fixture$seg, fixture$w)
    })
    expect_equal(result$data$kl_raw, vapply(raw, `[[`, numeric(1), "kl"))
    expect_equal(result$data$coverage_raw, vapply(raw, `[[`, numeric(1), "coverage"))
    rank <- order(if (priority == "higher") -x$distance_to_black else x$distance_to_black,
                   seq_len(nrow(x)))
    ordered <- optimize_sampling_distance(x[rank, ], fixture$seg, fixture$w, distances)
    expect_equal(result, ordered)
  }
  expect_error(optimize_sampling_distance(x, fixture$seg, fixture$w, distances,
    col_id = "missing", prioritize = "higher"), "missing")
})


test_that("HEALPix and Beckers propagate a fixed masked area reference", {
  geometry <- list(model = "radial_projection", parameters = equidistant_coef,
    center_of_vision = c(40, 30), normalization_radius = 30,
    normalization_angle = 90)
  g <- set_internal_geometry(c(80, 60), geometry)
  d <- direction_image(g, elevation = 32.5)
  w <- area_weight_image(g, geometry)
  m <- !is.na(d$zenith) & d$zenith <= 70
  for (seg in list(beckers_segmentation(d, 20), healpix_segmentation(d, 2))) {
    indices <- expand.grid(col = seq(1, 80, 3), row = seq(1, 60, 3))
    indices <- rem_unselected(indices, m & seg > 0)
    indices <- add_directions(indices, d = d)
    distances <- c(5, 20, 40)
    result <- suppressWarnings(optimize_sampling_distance(indices, seg, w, distances, m = m))
    expected <- lapply(distances, function(x) {
      retained <- rem_nearby_indices(indices, x, space = "spherical")
      assess_sampling_uniformity(retained, seg, w, m)
    })
    expect_equal(result$data$kl_raw, vapply(expected, `[[`, numeric(1), "kl"))
    expect_equal(result$data$coverage_raw, vapply(expected, `[[`, numeric(1), "coverage"))
  }
})
