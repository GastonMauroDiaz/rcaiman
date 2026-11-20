test_that("hull masks use image coordinates and ignore template values", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 6, ncols = 8))
  samples <- data.frame(x = c(1, 5, 5, 1), y_down = c(1, 1, 3, 3))
  m <- convex_hull_mask(samples, r)
  expected <- matrix(FALSE, nrow = 6, ncol = 8)
  expected[2:3, 2:5] <- TRUE
  expect_identical(terra::values(m, mat = FALSE), as.vector(t(expected)))
  expect_true(terra::compareGeom(m, r))
  expect_identical(names(m), "convex_hull")
  terra::values(r) <- NA_real_
  expect_identical(terra::values(convex_hull_mask(samples, r)), terra::values(m))
  # Duplicates, order, values and interior samples do not affect the hull.
  other <- rbind(samples[4:1, ], samples, data.frame(x = 2, y_down = 2))
  other$dn <- NA_real_
  expect_identical(terra::values(convex_hull_mask(other, r)), terra::values(m))
})

test_that("a nonrectangular hull selects cell centers without snapping samples", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 5, ncols = 5))
  samples <- data.frame(x = c(1.1, 4.1, 1.1), y_down = c(1.1, 1.1, 4.1))
  xy <- terra::xyFromCell(r, seq_len(ncell(r)))
  x <- xy[, 1]
  y_down <- nrow(r) - xy[, 2]
  expected <- x >= 1.1 & y_down >= 1.1 & x + y_down <= 5.2
  expect_identical(terra::values(convex_hull_mask(samples, r), mat = FALSE),
                   expected)
})

test_that("hulls outside the image are clipped rather than discarding samples", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 5, ncols = 5))
  samples <- data.frame(x = c(-10, 10, 10, -10),
                         y_down = c(-10, -10, 10, 10))
  expect_true(all(terra::values(convex_hull_mask(samples, r))))
  samples$x <- samples$x + 30
  expect_false(any(terra::values(convex_hull_mask(samples, r))))
})

test_that("hull masking composes with prediction and preserves interior values", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 8, ncols = 10))
  samples <- expand.grid(x = c(2, 4, 8), y_down = c(1, 3, 6))
  samples$dn <- samples$x - samples$y_down
  parameters <- fit_trend_surface(samples, np = 1)
  surface <- trend_surface_image(r, parameters)
  before <- terra::values(surface, mat = FALSE)
  m <- convex_hull_mask(samples, r)
  inside <- terra::values(m, mat = FALSE)
  surface[!m] <- NA
  observed <- terra::values(surface, mat = FALSE)
  expect_identical(observed[inside], before[inside])
  expect_true(all(is.na(observed[!inside])))
  expect_true(any(inside) && any(!inside))
  expect_equal(terra::values(trend_surface_image(r, parameters, m = m)),
               terra::values(surface))
})

test_that("degenerate hulls and undefined coordinates fail clearly", {
  r <- set_rcaiman_geometry(terra::rast(nrows = 5, ncols = 5))
  for (samples in list(data.frame(x = numeric(), y_down = numeric()),
                       data.frame(x = c(1, 2), y_down = c(1, 2)),
                       data.frame(x = c(1, 1, 1), y_down = c(1, 1, 1)),
                       data.frame(x = c(1, 2, 3), y_down = c(1, 2, 3)))) {
    expect_error(convex_hull_mask(samples, r), "non-collinear")
  }
  expect_error(convex_hull_mask(data.frame(x = c(1, 2, NA),
                                           y_down = c(1, 3, 4)), r), "NA")
  expect_error(convex_hull_mask(data.frame(col = 1:3, row = 1:3), r), "y_down")
})
