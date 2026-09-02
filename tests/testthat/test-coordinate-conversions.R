test_that("y and y_down use opposite vertical directions", {
  image_height <- 3
  y <- c(0, 0.5, 3)

  expect_equal(y_down_from_y(y, image_height), c(3, 2.5, 0))
  expect_equal(
    y_from_y_down(y_down_from_y(y, image_height), image_height),
    y
  )
})


test_that("coord_from_indices() returns cell centers", {
  indices <- data.frame(
    col = c(1, 2, 5),
    row = c(1, 3, 3)
  )

  expect_equal(
    coord_from_indices(indices$col, indices$row),
    data.frame(
      x = c(0.5, 1.5, 4.5),
      y_down = c(0.5, 2.5, 2.5)
    )
  )
})


test_that("indices_from_coord() assigns locations to containing cells", {
  coord <- data.frame(
    x = c(0, 0.5, 0.999, 1, 4.999),
    y_down = c(0, 0.5, 0.999, 1, 2.999)
  )

  expect_equal(
    indices_from_coord(coord$x, coord$y_down),
    data.frame(
      col = c(1L, 1L, 1L, 2L, 5L),
      row = c(1L, 1L, 1L, 2L, 3L)
    )
  )
})


test_that("xy_from_indices() returns cell centers with y increasing upward", {
  indices <- data.frame(
    col = c(1, 2, 5),
    row = c(1, 3, 3)
  )

  expect_equal(
    xy_from_indices(indices$col, indices$row, image_height = 3),
    data.frame(
      x = c(0.5, 1.5, 4.5),
      y = c(2.5, 0.5, 0.5)
    )
  )
})


test_that("cell centers round-trip between coord and indices", {
  indices <- data.frame(
    col = rep(1:5, times = 3),
    row = rep(1:3, each = 5)
  )
  coord <- coord_from_indices(indices$col, indices$row)

  expect_equal(
    indices_from_coord(coord$x, coord$y_down),
    indices
  )
})
