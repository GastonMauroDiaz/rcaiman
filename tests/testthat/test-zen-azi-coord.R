direction_images <- function(diameter = 10, orientation = 0) {
  z <- zenith_image(diameter, equidistant_coef)
  list(z = z, a = azimuth_image(z, orientation))
}


test_that("indices_from_zen_azi() characterizes direct cell assignment", {
  images <- direction_images()

  indices <- indices_from_zen_azi(
    images$z,
    images$a,
    zen = c(0, pi / 6, pi / 3),
    azi = c(0, pi / 2, pi)
  )

  expect_equal(
    indices,
    data.frame(
      col = c(6, 4, 6),
      row = c(6, 6, 9)
    )
  )
})


test_that("zen_azi_from_indices() characterizes cell lookup", {
  images <- direction_images()

  angles <- zen_azi_from_indices(
    images$z,
    images$a,
    col = c(5, 6, 8),
    row = c(5, 5, 6)
  )

  expect_equal(
    angles$zen,
    c(0.3141592654, 0.3141592654, 0.9424777961),
    tolerance = 1e-8
  )
  expect_equal(
    angles$azi,
    c(0.7853981634, 5.4977871438, 4.5149934205),
    tolerance = 1e-8
  )
})


test_that("indices_from_zen_azi() preserves duplicated cells", {
  images <- direction_images()

  indices <- indices_from_zen_azi(
    images$z,
    images$a,
    zen = c(pi / 6, pi / 6),
    azi = c(0, 0)
  )

  expect_equal(
    indices,
    data.frame(col = c(6, 6), row = c(4, 4))
  )
})


test_that("indices_from_zen_azi() preserves input order", {
  images <- direction_images()
  zen <- c(pi / 3, pi / 6, pi / 4)
  azi <- c(pi, pi / 2, 0)

  indices <- indices_from_zen_azi(images$z, images$a, zen, azi)
  reversed_indices <- indices_from_zen_azi(
    images$z,
    images$a,
    rev(zen),
    rev(azi)
  )

  expect_equal(
    indices,
    data.frame(
      col = c(6, 4, 6),
      row = c(9, 6, 3)
    )
  )
  expect_identical(
    reversed_indices,
    data.frame(
      col = c(6, 4, 6),
      row = c(3, 6, 9)
    )
  )
})


test_that("indices_from_zen_azi() returns NA outside the raster", {
  images <- direction_images()

  indices <- indices_from_zen_azi(
    images$z,
    images$a,
    zen = pi,
    azi = 0
  )

  expect_equal(indices, data.frame(col = NA_integer_, row = NA_integer_))
})


test_that("indices_from_zen_azi() preserves mixed angular correspondence", {
  images <- direction_images()
  zen <- c(pi, pi / 6)
  azi <- c(0, pi / 2)

  observed <- indices_from_zen_azi(images$z, images$a, zen, azi)
  expected <- do.call(
    rbind,
    lapply(seq_along(zen), function(i) {
      indices_from_zen_azi(images$z, images$a, zen[i], azi[i])
    })
  )

  expect_equal(observed, expected)
})


test_that("zen_azi_from_indices() handles image borders", {
  images <- direction_images()

  angles <- zen_azi_from_indices(
    images$z,
    images$a,
    col = c(1, 10),
    row = c(1, 10)
  )

  expect_true(all(is.na(angles)))
  expect_error(
    zen_azi_from_indices(images$z, images$a, 11, 1),
    "`col` must be up to 10",
    fixed = TRUE
  )
  expect_error(
    zen_azi_from_indices(images$z, images$a, 1, 11),
    "`row` must be up to 10",
    fixed = TRUE
  )
  expect_error(zen_azi_from_indices(images$z, images$a, 0, 1))
  expect_error(zen_azi_from_indices(images$z, images$a, 1, 0))
})


test_that("angular conversion validates vector lengths and ranges", {
  images <- direction_images()

  expect_error(
    zen_azi_from_indices(images$z, images$a, c(1, 2), 1),
    "must have same length"
  )
  expect_error(
    indices_from_zen_azi(images$z, images$a, c(0, pi / 2), 0),
    "must have same length"
  )
  expect_error(indices_from_zen_azi(images$z, images$a, -0.1, 0))
  expect_error(indices_from_zen_azi(images$z, images$a, pi + 0.1, 0))
  expect_error(indices_from_zen_azi(images$z, images$a, 0, 2 * pi + 0.1))
})


test_that("angular conversion honors azimuth image orientation", {
  skip_if_not_installed("imager")
  images <- direction_images(100, orientation = 45)
  indices <- data.frame(col = 51, row = 26)

  angles <- zen_azi_from_indices(
    images$z,
    images$a,
    indices$col,
    indices$row
  )
  recovered_indices <- indices_from_zen_azi(
    images$z,
    images$a,
    angles$zen,
    angles$azi
  )

  expect_equal(recovered_indices, indices)
})
