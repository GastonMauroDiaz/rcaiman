geometry_test_image_path <- function() {
  system.file("external/DSCN4500.JPG", package = "rcaiman")
}


test_that("read_caim() and crop_caim() compose on a real image", {
  path <- geometry_test_image_path()
  full_caim <- read_caim(path, fill_na = "none")

  observed <- crop_caim(
    full_caim,
    snapped_zenith_coord = c(5, 5),
    horizon_radius = 3
  )
  expected <- crop_caim(
    full_caim,
    upper_left_indices = c(3, 3),
    width = 6,
    height = 6
  )

  expect_equal(dim(observed), c(6, 6, 3))
  expect_equal(terra::values(observed), terra::values(expected))
})


test_that("expand_non_circular() preserves cells around the snapped zenith", {
  caim <- crop_caim(
    read_caim(geometry_test_image_path()),
    upper_left_indices = c(1000, 800),
    width = 6,
    height = 4
  )
  z <- zenith_image(8, equidistant_coef)

  observed <- expand_non_circular(
    caim,
    z,
    snapped_zenith_coord = c(3, 2)
  )
  target_cells <- terra::cellFromRowCol(
    z,
    row = rep(3:6, each = 6),
    col = rep(2:7, times = 4)
  )

  expect_equal(dim(observed), c(8, 8, 3))
  expect_equal(
    unname(terra::values(observed)[target_cells, ]),
    unname(terra::values(caim))
  )
  expect_true(all(is.na(terra::values(observed)[1, ])))
})


test_that("expand_non_circular() uses z when the translated image is wider", {
  caim <- set_rcaiman_geometry(terra::rast(nrows = 4, ncols = 10))
  terra::values(caim) <- seq_len(terra::ncell(caim))
  names(caim) <- "dn"
  z <- zenith_image(8, equidistant_coef)

  observed <- expand_non_circular(
    caim,
    z,
    snapped_zenith_coord = c(5, 2)
  )
  source_matrix <- matrix(
    terra::values(caim)[, 1],
    nrow = 4,
    byrow = TRUE
  )
  expected_matrix <- matrix(NA_real_, nrow = 8, ncol = 8)
  expected_matrix[3:6, ] <- source_matrix[, 2:9]

  expect_true(terra::compareGeom(observed, z))
  expect_equal(dim(observed), c(8, 8, 1))
  expect_equal(
    terra::values(observed)[, 1],
    as.vector(t(expected_matrix))
  )
})


test_that("expand_non_circular() requires a snapped grid intersection", {
  caim <- crop_caim(
    read_caim(geometry_test_image_path()),
    upper_left_indices = c(1000, 800),
    width = 6,
    height = 4
  )
  z <- zenith_image(8, equidistant_coef)

  expect_error(
    expand_non_circular(caim, z, snapped_zenith_coord = c(3.5, 2)),
    "integer-like"
  )
})


test_that("expanded points identify the same cells as the expanded image", {
  caim <- crop_caim(
    read_caim(geometry_test_image_path()),
    upper_left_indices = c(1000, 800),
    width = 6,
    height = 4
  )
  z <- zenith_image(8, equidistant_coef)
  snapped_zenith_coord <- c(3, 2)
  original_indices <- data.frame(
    col = c(1, 3, 6),
    row = c(1, 2, 4)
  )
  original_coord <- coord_from_indices(
    original_indices$col,
    original_indices$row
  )

  expanded_coord <- expand_non_circular_points(
    original_coord,
    z,
    snapped_zenith_coord
  )
  expanded_indices <- indices_from_coord(
    expanded_coord$x,
    expanded_coord$y_down
  )
  expanded_caim <- expand_non_circular(
    caim,
    z,
    snapped_zenith_coord
  )
  original_cells <- terra::cellFromRowCol(
    caim,
    original_indices$row,
    original_indices$col
  )
  expanded_cells <- terra::cellFromRowCol(
    z,
    expanded_indices$row,
    expanded_indices$col
  )

  expect_equal(expanded_coord$x, original_coord$x + 1)
  expect_equal(expanded_coord$y_down, original_coord$y_down + 2)
  expect_equal(
    unname(terra::values(expanded_caim)[expanded_cells, ]),
    unname(terra::values(caim)[original_cells, ])
  )
})


test_that("non-circular point expansion is exactly reversible", {
  z <- zenith_image(8, equidistant_coef)
  original_coord <- data.frame(
    x = c(0.5, 0.5, 7.5),
    y_down = c(0.5, 0.5, 5.5),
    id = c("duplicate-a", "duplicate-b", "outside")
  )

  expanded_coord <- expand_non_circular_points(
    original_coord,
    z,
    snapped_zenith_coord = c(3, 2)
  )
  expanded_coord <- add_indices(expanded_coord)
  expanded_coord <- add_xy(expanded_coord, image_height = nrow(z))
  expanded_coord$relative_radius <- c(0.2, 0.2, 0.9)
  expanded_coord$image_angle <- c(0, 0, pi / 2)
  observed <- unexpand_non_circular_points(
    expanded_coord,
    z,
    snapped_zenith_coord = c(3, 2)
  )

  expect_identical(observed, original_coord)
  expect_identical(expanded_coord$id, original_coord$id)
})


test_that("non-circular point expansion removes previous-frame geometry", {
  z <- zenith_image(8, equidistant_coef)
  original_coord <- data.frame(
    x = c(0.5, 2.5),
    y_down = c(0.5, 1.5),
    col = c(1L, 3L),
    row = c(1L, 2L),
    y = c(5.5, 4.5),
    relative_radius = c(0.2, 0.4),
    image_angle = c(0, pi / 2),
    zen = c(0.3, 0.6),
    azi = c(0.4, 0.8),
    id = c("a", "b")
  )

  observed <- expand_non_circular_points(
    original_coord,
    z,
    snapped_zenith_coord = c(3, 2)
  )

  expect_identical(
    names(observed),
    c("x", "y_down", "zen", "azi", "id")
  )
  expect_equal(observed$x, original_coord$x + 1)
  expect_equal(observed$y_down, original_coord$y_down + 2)
  expect_identical(observed$zen, original_coord$zen)
  expect_identical(observed$azi, original_coord$azi)
  expect_identical(observed$id, original_coord$id)
})


test_that("non-circular point expansion accepts a structured empty input", {
  z <- zenith_image(8, equidistant_coef)
  original_coord <- data.frame(
    x = numeric(),
    y_down = numeric(),
    id = character()
  )

  expanded_coord <- expand_non_circular_points(
    original_coord,
    z,
    snapped_zenith_coord = c(3, 2)
  )

  expect_identical(expanded_coord, original_coord)
  expect_identical(
    unexpand_non_circular_points(
      expanded_coord,
      z,
      snapped_zenith_coord = c(3, 2)
    ),
    original_coord
  )
})
