test_that("mask_image() works", {
  local_edition(3)
  skip_on_cran()

  z <- zenith_image(1000, lens())
  a <- azimuth_image(z)

  write_img_and_return_path <- function() {
    r <- mask_image(z, a, c(20, 70), c(90, 180))
    path <- tempfile(fileext = ".tif")
    write_bin(r, path)
    path
  }
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "mask_zlim_alim.tif")
  file.remove(path)

  write_img_and_return_path <- function() {
    r <- mask_image(z, zlim = c(20, 70))
    path <- tempfile(fileext = ".tif")
    write_bin(r, path)
    path
  }
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "mask_zlim.tif")
  file.remove(path)

  write_img_and_return_path <- function() {
    r <- mask_image(z)
    path <- tempfile(fileext = ".tif")
    write_bin(r, path)
    path
  }
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "mask_default.tif")
  file.remove(path)

})

test_that("masking() works", {
  local_edition(3)
  skip_on_cran()

  r <- read_caim()
  z <- zenith_image(ncol(r), lens())
  a <- azimuth_image(z)

  write_img_and_return_path <- function() {
    m <- mask_image(z, a, c(20, 70), c(90, 180))
    r <- normalize(r, 0, 255) %>% masking(., m)
    path <- tempfile(fileext = ".tif")
    write_caim(r * 2^8, path, 8)
    path
  }
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "masking_caim.tif")
  file.remove(path)

  write_img_and_return_path <- function() {
    m <- mask_image(z, a, c(20, 70), c(90, 180))
    r <- apply_thr(r$Blue, 125)
    r <- masking(r, m)
    path <- tempfile(fileext = ".tif")
    write_caim(r * 2^8, path, 8)
    path
  }
  path <- write_img_and_return_path()
  expect_snapshot_file(path, "masking_bin.tif")
  file.remove(path)

})


