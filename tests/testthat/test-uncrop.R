test_that("image and internal geometry return to the same original cells", {
  image <- set_rcaiman_geometry(terra::rast(ncols = 8, nrows = 6, nlyrs = 3))
  values <- cbind(1:48, 101:148, 201:248)
  values[19, 2] <- NA_real_
  terra::values(image) <- values
  g <- set_internal_geometry(c(8, 6), list(
    center_of_vision = c(4, 3),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = 3,
    normalization_angle = 90
  ))
  original_dim <- c(ncol(image), nrow(image))
  original_image <- terra::values(image)
  original_g <- terra::values(g)

  # Interior, corners, a single cell, and the entire grid.
  frames <- list(c(2, 2, 4, 3), c(1, 1, 3, 2), c(6, 1, 3, 2),
                 c(1, 5, 3, 2), c(6, 5, 3, 2), c(8, 6, 1, 1),
                 c(1, 1, 8, 6))
  for (frame in frames) {
    upper_left <- frame[1:2]
    width <- frame[3]
    height <- frame[4]
    image_crop <- crop_caim(image, upper_left, width, height)
    g_crop <- crop_caim(g, upper_left, width, height)
    image_back <- uncrop_caim(image_crop, upper_left, original_dim)
    g_back <- uncrop_caim(g_crop, upper_left, original_dim)
    expect_true(terra::compareGeom(image_crop, g_crop))
    expect_true(terra::compareGeom(image_back, g_back, image))

    cells <- terra::cellFromRowCol(
      image,
      row = rep(seq(upper_left[2], length.out = height), each = width),
      col = rep(seq(upper_left[1], length.out = width), times = height)
    )
    originals <- list(image, g)
    crops <- list(image_crop, g_crop)
    restored <- list(image_back, g_back)
    for (i in seq_along(originals)) {
      expected <- terra::values(originals[[i]])
      expected[-cells, ] <- NA_real_
      expect_equal(terra::values(restored[[i]]), expected)
      expect_identical(names(restored[[i]]), names(originals[[i]]))
      expect_equal(unname(as.vector(terra::ext(crops[[i]]))),
                   c(0, width, 0, height))
    }
  }
  expect_equal(terra::values(image), original_image)
  expect_equal(terra::values(g), original_g)
  expect_equal(unname(as.vector(terra::ext(image))), c(0, 8, 0, 6))
})

test_that("uncropping angular geometry preserves external directions", {
  g <- set_internal_geometry(c(8, 6), list(
    center_of_vision = c(4, 3),
    model = "radial_projection",
    parameters = 2 / pi,
    normalization_radius = 3,
    normalization_angle = 90
  ))
  cropped <- crop_caim(g, c(2, 2), 4, 3)
  before <- uncrop_caim(direction_image(cropped, 90, 45, 32.5),
                        c(2, 2), c(8, 6))
  after <- direction_image(uncrop_caim(cropped, c(2, 2), c(8, 6)),
                            90, 45, 32.5)
  expect_equal(terra::values(before), terra::values(after))
  expect_identical(names(before), names(after))
})

test_that("uncrop_caim() handles a single layer and an all-NA crop", {
  r <- set_rcaiman_geometry(terra::rast(ncols = 2, nrows = 3))
  terra::values(r) <- NA_real_
  names(r) <- "dn"
  observed <- uncrop_caim(r, c(3, 2), c(5, 6))
  expect_equal(dim(observed), c(6, 5, 1))
  expect_identical(names(observed), "dn")
  expect_true(all(is.na(terra::values(observed))))
})

test_that("uncrop_caim() checks placement and original grid dimensions", {
  r <- set_rcaiman_geometry(terra::rast(ncols = 2, nrows = 3))
  terra::values(r) <- 1:6
  expect_error(uncrop_caim(r, c(4, 1), c(4, 4)), "extends outside")
  expect_error(uncrop_caim(r, c(1, 3), c(4, 4)), "extends outside")
  expect_error(uncrop_caim(r, c(0, 1), c(4, 4)), "positive")
  expect_error(uncrop_caim(r, c(1.5, 1), c(4, 4)), "integer-like")
  expect_error(uncrop_caim(r, c(1, 1), c(4, 4.5)), "integer-like")
  expect_error(uncrop_caim(r, c(1, 1), c(4, 0)), "positive")
  expect_error(uncrop_caim(r, c(1, 1), c(4, NA)), "NA")
  expect_error(uncrop_caim(r, c(1, 1), c(4, Inf)), "integer-like")
  expect_error(uncrop_caim(r, c(1, 1), NULL), "NULL")
})
