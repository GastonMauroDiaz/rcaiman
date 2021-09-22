# \dontrun{
# path <- system.file("external", package = "rcaiman")
# my_file <- paste0(path, "/DSCN5548.JPG")
# download.file("https://osf.io/kp7rx/download", my_file,
#                method = "auto", mode = "wb"
# )
#
# r <- read_caim(file.path(path, "DSCN5548.JPG"),
#                c(1280, 960) - 745, 745 * 2, 745 * 2)
# z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
# a <- azimuth_image(z)
# thr <- autothresholdr::auto_thresh(r$Blue[], "IsoData")
# bin <- apply_thr(r$Blue, thr[1] * 1.25)
# blue <- gbc(r$Blue)
# sky <- model_sky_dn(blue, z, a, bin, parallel = FALSE)
# path <- tempfile(fileext = ".tif")
# write_caim(sky$image * 2^8, path, 8)
# }
