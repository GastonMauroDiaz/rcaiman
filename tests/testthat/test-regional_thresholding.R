test_that("regional_thresholding() works", {
  local_edition(3)
  skip_on_cran()
  write_regional_bin_and_return_path <- function(r) {
    r <- read_caim()
    blue <- gbc(r$Blue)
    z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
    rings <- rings_segmentation(z, 10)
    bin <- regional_thresholding(blue, rings, "IsoData")
    path <- tempfile(fileext = ".tif")
    write_bin(bin, path)
    path
  }
  path <- write_regional_bin_and_return_path()
  expect_snapshot_file(path, "isodata_regional.tif")
  file.remove(path)

  write_regional_bin_and_return_path <- function(r) {
    r <- read_caim()
    blue <- gbc(r$Blue)
    z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
    rings <- rings_segmentation(z, 10)
    bin <- regional_thresholding(blue, rings, "mblt_default")
    path <- tempfile(fileext = ".tif")
    write_bin(bin, path)
    path
  }
  path <- write_regional_bin_and_return_path()
  expect_snapshot_file(path, "mblt_regional.tif")
  file.remove(path)
})
