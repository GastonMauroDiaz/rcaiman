#' Build a sky image
#'
#' Build an above canopy image from a single below canopy image.
#'
#' @inheritParams ootb_mblt
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
build_sky <- function(r, z, a, complicado = FALSE) {
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop(paste("Package \"lidR\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }
  model <- autofit_cie_sky_model(r, z, a)

  sky_cie <- model$relative_luminance * model$zenith_dn
  residu <- r - sky_cie

  bin <- model$mblt$bin

  ##############
  if (!complicado) {

  rcells <- r
  rcells[] <- 1:ncell(r)
  xy <- xyFromCell(r, rcells[bin])
  las <- .make_fake_las(xy[, 1], xy[, 2], residu[bin])
  las@data$Classification <- 2
  residu <- lidR::grid_terrain(las, res = 1,
                               algorithm = lidR::knnidw())
  residu <- resample(residu, r)

  }
  ###############
  if (complicado) {

  m1 <- mask_image(z, a, alim = c(0,30)) & bin
  m2 <- mask_image(z, a, alim = c(330,360)) & bin
  m3 <- mask_image(z, zlim = c(0,30)) & bin
  las <- .make_fake_las(c(a[bin]     , a[m1] + 360, a[m2] - 360, abs(a[m3] - 360)),
                        c(z[bin]     , z[m1]      , z[m2]      , abs(z[m3] - 30) - 30),
                        c(residu[bin], residu[m1] , residu[m2] , residu[m3])
                        )
  las@data$Classification <- 2

  suppressWarnings(
    suppressMessages(residu <- lidR::grid_terrain(
                                            las, res = 1,
                                            # algorithm = lidR::tin()
                                            algorithm = lidR::knnidw(k = 10,
                                                                     p = 1,
                                                                     rmax = 20)
                                            )
                   ))
  residu <- crop(residu, extent(0,360,0,90))
  g <- sky_grid_segmentation(z, a, 1)
  az <- xyFromCell(residu, 1:ncell(residu)) + 0.5
  g_id <- az[,1] * 1000 + az[,2]
  rcl <- data.frame(g_id, residu[])
  residu <- reclassify(g, rcl)
  residu[residu > 1000] <- NA

  }
  #############################

  sky <- sky_cie + residu
  sky <- cover(sky, sky_cie)

  list(model = model,
       interpolated_residual = residu,
       sky = sky)
}
