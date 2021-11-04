#' Build a sky image
#'
#' Build an above canopy image from a single below canopy image.
#'
#' @inheritParams ootb_mblt
#' @inheritParams lidR::knnidw
#' @return
#' @export
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' sky <- build_sky(blue, z, a)
#' plot(sky$sky)
#' }
build_sky <- function(r, z, a, k = 50, p = 1, rmax = 200) {
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop(paste("Package \"lidR\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }
  model <- autofit_cie_sky_model(r, z, a)

  sky_cie <- model$relative_luminance * model$zenith_dn
  bin <- suppressWarnings(apply_thr(r, thr_image(sky_cie, 0, 0.5)))

  flat <- residu <- sky_cie - r
  flat[] <- 0
  residu <- focal(residu, matrix(1, ncol = 3, nrow = 3), fun = max)
  m <- mask_image(z)
  residu_s <- suppressWarnings(fit_trend_surface(residu, m, bin, flat,
                                                 fact = 1))
  sky_s <- sky_cie - residu_s$image

  rcells <- r
  rcells[] <- 1:ncell(r)
  xy <- xyFromCell(r, rcells[bin])
  las <- .make_fake_las(xy[, 1], xy[, 2], residu[bin])
  las@data$Classification <- 2

  residu <- lidR::grid_terrain(las, res = 1,
                                       algorithm = lidR::knnidw(k, p, rmax))

  residu <- resample(residu, r)
  sky <- sky_cie - residu
  sky <- cover(sky, sky_cie)
  list(model = model,
       interpolated_residual = residu,
       sky = sky)
}
