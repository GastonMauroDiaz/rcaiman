ootb_sky_reconstruction <- function(r, z, a, lens_coef, is_open_forest = FALSE) {
  g <- sky_grid_segmentation(z, a, 10)
  bin <- ootb_mblt(r, z, a, is_open_forest)$bin
  if (!is_open_forest) {
    m <- mask_image(z, zlim = c(70, 90))
    bin[m] <- 0
  }
  sky_marks <- extract_sky_marks(r, bin, g)
  sun_mark <- extract_sun_mark(r, bin, z, a, g)
  model <- fit_cie_sky_model(r, z, a, sky_marks, sun_mark)
  sky_cie <- model$relative_luminance * model$zenith_dn
  residu <- sky_cie - r
  residu_i <- interpolate_dns(residu, z, a, sky_marks, lens_coef)
  sky <- sky_cie - residu_i
  sky <- cover(sky, sky_cie)
}
