test_that("read_points() reads coord and preserves additional columns", {
  path <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      x = c(0.5, 12.75),
      y_down = c(1.5, 8.25),
      point_id = c("a", "b")
    ),
    path,
    row.names = FALSE
  )

  observed <- read_points(path)

  expect_identical(names(observed), c("x", "y_down", "point_id"))
  expect_equal(observed$x, c(0.5, 12.75))
  expect_equal(observed$y_down, c(1.5, 8.25))
  expect_identical(observed$point_id, c("a", "b"))
})

test_that("read_points() rejects cell indices", {
  path <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(col = c(1, 2), row = c(1, 2)),
    path,
    row.names = FALSE
  )

  expect_error(
    read_points(path),
    "required columns: x, y_down",
    fixed = TRUE
  )
})

test_that("write_points() writes coord first and preserves other columns", {
  path <- tempfile(fileext = ".csv")
  sampling_coord <- data.frame(
    point_id = c("a", "b"),
    y_down = c(1.5, 8.25),
    value = c(10, 20),
    x = c(0.5, 12.75)
  )

  write_points(sampling_coord = sampling_coord, path = path)
  observed <- read.csv(path)

  expect_identical(
    names(observed),
    c("x", "y_down", "point_id", "value")
  )
  expect_equal(observed$x, sampling_coord$x)
  expect_equal(observed$y_down, sampling_coord$y_down)
  expect_identical(observed$point_id, sampling_coord$point_id)
  expect_equal(observed$value, sampling_coord$value)
})

test_that("write_points() rejects cell indices", {
  path <- tempfile(fileext = ".csv")

  expect_error(
    write_points(data.frame(col = 1, row = 1), path),
    "required columns: x, y_down",
    fixed = TRUE
  )
})
