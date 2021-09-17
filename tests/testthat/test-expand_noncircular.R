test_that("expand_noncircular() works", {
  local_edition(3)
  skip_on_cran()

  read_and_expand_DSC_2881 <- function() {
    path <- system.file("external", package = "rcaiman")
    my_file <- paste0(path, "/DSC_2881.JPG")

    if (!file.exists(my_file)) {
      download.file("https://osf.io/x8urg/download", my_file,
                    method = "auto", mode = "wb"
      )
    }

    r <- read_caim(file.path(path, "DSC_2881.JPG"))

    diameter <- calc_diameter(lens("Nikkor_10.5_mm"), 1202, 53)
    zenith_colrow <- c(1503, 998)
    z <- zenith_image(diameter, lens("Nikkor_10.5_mm"))
    r <- expand_noncircular(r, z, zenith_colrow)

    path <- tempfile(fileext = ".tif")
    writeRaster(r, path, datatype = "INT1U", overwrite = TRUE)

    path
  }


  path <- read_and_expand_DSC_2881()
  expect_snapshot_file(path, "DSC_2881_expanded.tif")
  file.remove(path)
})
