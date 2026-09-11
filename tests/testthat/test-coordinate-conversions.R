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


test_that("add_coord() appends cell-center coordinates", {
  sampling_indices <- data.frame(
    col = c(1L, 5L, 1L),
    row = c(1L, 3L, 1L),
    id = c("first", "second", "duplicate"),
    row.names = c("a", "b", "c")
  )

  observed <- add_coord(sampling_indices)

  expect_identical(
    names(observed),
    c("col", "row", "id", "x", "y_down")
  )
  expect_equal(observed$x, c(0.5, 4.5, 0.5))
  expect_equal(observed$y_down, c(0.5, 2.5, 0.5))
  expect_identical(observed$id, sampling_indices$id)
  expect_identical(row.names(observed), row.names(sampling_indices))
})


test_that("add_indices() appends containing-cell indices", {
  sampling_coord <- data.frame(
    x = c(0, 0.999, 1, 4.999),
    y_down = c(0.5, 0.999, 1, 2.999),
    id = letters[1:4]
  )

  observed <- add_indices(sampling_coord)

  expect_identical(
    names(observed),
    c("x", "y_down", "id", "col", "row")
  )
  expect_identical(observed$col, c(1L, 1L, 2L, 5L))
  expect_identical(observed$row, c(1L, 1L, 2L, 3L))
  expect_identical(observed$id, sampling_coord$id)
})


test_that("add_xy() appends the upward vertical coordinate", {
  sampling_coord <- data.frame(
    x = c(0, 0.5, 3),
    y_down = c(0, 0.5, 3),
    id = letters[1:3]
  )

  observed <- add_xy(sampling_coord, image_height = 3)

  expect_identical(names(observed), c("x", "y_down", "id", "y"))
  expect_equal(observed$y, c(3, 2.5, 0))
  expect_identical(observed[names(sampling_coord)], sampling_coord)
})


test_that("enrichment helpers preserve structured zero-row inputs", {
  empty_indices <- data.frame(
    col = integer(),
    row = integer(),
    id = character()
  )
  empty_coord <- data.frame(
    x = numeric(),
    y_down = numeric(),
    id = character()
  )

  expect_identical(
    names(add_coord(empty_indices)),
    c("col", "row", "id", "x", "y_down")
  )
  expect_identical(
    names(add_indices(empty_coord)),
    c("x", "y_down", "id", "col", "row")
  )
  expect_identical(
    names(add_xy(empty_coord, image_height = 3)),
    c("x", "y_down", "id", "y")
  )
})


test_that("enrichment helpers reject missing and existing target columns", {
  expect_error(
    add_coord(data.frame(col = 1L)),
    "required columns: row",
    fixed = TRUE
  )
  expect_error(
    add_coord(data.frame(col = 1L, row = 1L, x = 0.5)),
    "reserved columns: x",
    fixed = TRUE
  )
  expect_error(
    add_indices(data.frame(x = 0.5, y_down = 0.5, row = 1L)),
    "reserved columns: row",
    fixed = TRUE
  )
  expect_error(
    add_xy(data.frame(x = 0.5, y_down = 0.5, y = 2.5), 3),
    "reserved columns: y",
    fixed = TRUE
  )
})


test_that("keep_coord() retains coord and frame-independent columns", {
  sampling_coord <- data.frame(
    id = c("first", "second"),
    x = c(0.5, 2.5),
    y_down = c(1.5, 3.5),
    col = c(1L, 3L),
    row = c(2L, 4L),
    y = c(4.5, 2.5),
    relative_radius = c(0.1, 0.5),
    image_angle = c(0, pi / 2),
    zen = c(0.2, 0.8),
    azi = c(0.4, 1.2),
    dn = c(10, 20)
  )

  observed <- keep_coord(sampling_coord)
  expected <- sampling_coord[c("x", "y_down", "id", "zen", "azi", "dn")]

  expect_identical(observed, expected)
  expect_identical(keep_coord(observed), observed)
})


test_that("keep_coord() preserves a structured zero-row coord", {
  sampling_coord <- data.frame(
    id = character(),
    x = numeric(),
    y_down = numeric(),
    col = integer(),
    row = integer(),
    y = numeric(),
    zen = numeric()
  )

  observed <- keep_coord(sampling_coord)

  expect_identical(names(observed), c("x", "y_down", "id", "zen"))
  expect_equal(nrow(observed), 0)
})
