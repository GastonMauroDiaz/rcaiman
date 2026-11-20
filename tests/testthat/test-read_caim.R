test_that("reading an RGB image and cropping preserve bands and cells", {
  path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
  r <- read_caim(path, fill_na = "none")
  expect_setequal(names(r), c("Red", "Green", "Blue"))
  cropped <- crop_caim(r, c(3, 3), width = 6, height = 6)
  cells <- terra::cellFromRowCol(r, rep(3:8, each = 6), rep(3:8, times = 6))
  expect_equal(dim(cropped), c(6, 6, 3))
  expect_equal(terra::values(cropped), terra::values(r)[cells, ])
})
