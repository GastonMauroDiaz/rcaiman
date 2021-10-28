test_that("fit_cone_shaped_model() works", {
  get_sky_from_DSCN5547 <- function() {
    path <- system.file("external", package = "rcaiman")
    my_file <- paste0(path, "/DSCN5548.JPG")

    if (!file.exists(my_file)) {
      download.file("https://osf.io/kp7rx/download", my_file,
                    method = "auto", mode = "wb"
      )
    }

    r <- read_caim(file.path(path, "DSCN5548.JPG"),
                   c(1280, 960) - 745, 745 * 2, 745 * 2)
    z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
    a <- azimuth_image(z)
    thr <- autothresholdr::auto_thresh(r$Blue[], "IsoData")
    bin <- apply_thr(r$Blue, thr[1] * 1.25)
    blue <- gbc(r$Blue)
    sky <- fit_cone_shaped_model(blue, z, a, bin, parallel = FALSE)
    path <- tempfile(fileext = ".tif")
    write_caim(sky$image * 2^8, path, 8)

    path
  }

  local_edition(3)
  skip_on_cran()
  path <- get_sky_from_DSCN5547()
  expect_snapshot_file(path, "sky_from_DSCN5547_only.tif")
  file.remove(path)
})
