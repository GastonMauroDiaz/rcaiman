test_that("calc_relative_radius returns between the 0-to-1 range", {
  expect_equal(calc_relative_radius(0, lens("equiangular")), 0)
  expect_equal(calc_relative_radius(90, lens("equiangular")), 1)
  expect_equal(calc_relative_radius(0, lens("Nikon_FCE9")), 0)
  expect_equal(calc_relative_radius(90, lens("Nikon_FCE9")), 1)
  expect_equal(calc_relative_radius(0, lens("Nikkor_10.5_mm")), 0)
  expect_equal(calc_relative_radius(90, lens("Nikkor_10.5_mm")), 1)
  expect_equal(calc_relative_radius(0, lens("Soligor_fisheye")), 0)
  expect_equal(calc_relative_radius(90, lens("Soligor_fisheye")), 1)
  expect_equal(calc_relative_radius(0, lens("Olloclip")), 0)
  expect_equal(calc_relative_radius(90, lens("Olloclip")), 1)
})

test_that("relative_radius_image works", {
  write_img_and_return_path <- function() {
    r <- relative_radius_image(1490)
    path <- tempfile(fileext = ".tif")
    suppressWarnings(writeRaster(r, path, overwrite = TRUE))
    path
  }
  local_edition(3)
  skip_on_cran()
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "relative_radius_image.tif")
  file.remove(path)
})

test_that("zenith_image works", {
  write_img_and_return_path <- function() {
    r <- zenith_image(1490, lens_coef = lens("Nikon_FCE9"))
    path <- tempfile(fileext = ".tif")
    suppressWarnings(writeRaster(r, path, overwrite = TRUE))
    path
  }
  local_edition(3)
  skip_on_cran()
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "zenith_image.tif")
  file.remove(path)
})

test_that("azimuth_image works", {
  write_img_and_return_path <- function() {
    z <- zenith_image(1490, lens_coef = lens("Nikon_FCE9"))
    r <- azimuth_image(z)
    path <- tempfile(fileext = ".tif")
    suppressWarnings(writeRaster(r, path, overwrite = TRUE))
    path
  }
  local_edition(3)
  skip_on_cran()
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "azimuth_image.tif")
  file.remove(path)
})
