

test_that("a call to read_caim() returs the example photo", {
  write_example_photo_and_return_path <- function(r) {
    r <- read_caim()
    path <- tempfile(fileext = ".tif")
    writeRaster(r, path, datatype = "INT1U", overwrite = TRUE)
    path
  }
  local_edition(3)
  skip_on_cran()
  path <- write_example_photo_and_return_path()
  expect_snapshot_file(path, "example.tif")
  file.remove(path)
})


test_that("a call to read_caim() returs a RasterBrick", {
  expect_is(r <- read_caim(), "RasterBrick")
  expect_equal(nlayers(r), 3)
})

test_that("read_caim() assigns good layer names", {
  expect_setequal(names(read_caim()), c("Red", "Green", "Blue"))
})




test_that("a ROI from a photo can be read correctly", {
  read_19_crop_save_and_return_path <- function() {
    path <- system.file("external", package = "rcaiman")
    my_file <- paste0(path, "/19.JPG")

    if (!file.exists(my_file)) {
      download.file("https://osf.io/zh5md/download", my_file,
        method = "auto", mode = "wb"
      )
    }

    r <- read_caim(file.path(path, "19.JPG"),
      upper_left = c(534, 237),
      width = 1530, height = 1476
    )
    path <- tempfile(fileext = ".tif")
    writeRaster(r, path, datatype = "INT1U", overwrite = TRUE)
    path
  }
  local_edition(3)
  skip_on_cran()
  path <- read_19_crop_save_and_return_path()
  expect_snapshot_file(path, "19_cropped.tif")
  file.remove(path)
})
