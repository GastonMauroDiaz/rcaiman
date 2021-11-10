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
interpolate_dns <- function(r, z, a, sky_marks, lens_coef, use_window = TRUE) {
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop(paste("Package \"lidR\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }

  cells <- cellFromRowCol(a, sky_marks$row, sky_marks$col)
  sky_marks$a <- a[cells]
  sky_marks$z <- z[cells]

  if (use_window) {
    xy <-  xyFromCell(r, cells)
    sky_marks$dn <-  extract(r, xy, buffer = 1.5, fun = "mean")
  } else {
    sky_marks$dn <-  r[cells]
  }

  # las <- .make_fake_las(xy[, 1], xy[, 2], sky_marks$dn)
  # las@data$Classification <- 2
  # ir <- lidR::grid_terrain(las, res = 1,
  #                              algorithm = lidR::knnidw())
  # ir <- resample(ir, r)

  m1 <- sky_marks$a > 0 & sky_marks$a < 30
  m2 <- sky_marks$a > 330 & sky_marks$a < 360
  m3 <- sky_marks$z > 0 & sky_marks$z < 30
  las <- .make_fake_las(c(sky_marks$a , sky_marks$a[m1] + 360, sky_marks$a[m2] - 360, abs(sky_marks$a[m3] - 360)),
                        c(sky_marks$z , sky_marks$z[m1]      , sky_marks$z[m2]      , abs(sky_marks$z[m3] - 30) - 30),
                        c(sky_marks$dn, sky_marks$dn[m1]     , sky_marks$dn[m2]     , sky_marks$dn[m3])
                        )
  las@data$Classification <- 2

  temp <- raster()
  projection(temp) <- NA
  extent(temp) <- extent(-30, 390, -30, 90)
  .res <- (ncell(temp) / ncell(r)) * 10
  suppressWarnings(
    suppressMessages(ir <- lidR::grid_terrain(
                                            las, res = .res,
                                            # algorithm = lidR::tin()
                                            algorithm = lidR::knnidw(k = 10,
                                                                     p = 1,
                                                                     rmax = 20)
                                            )
                   ))
  ir <- crop(ir, extent(0,360,0,90))
  az <- xyFromCell(ir, 1:ncell(ir))
  rr <- calc_relative_radius(az[,2], lens_coef)
  pol <- data.frame(theta = az[,1] * pi/180 + pi/2,
                    r = rr * 90 * pi/180,
                    z = ir[])
  cart <- pracma::pol2cart(as.matrix(pol))
  p <- SpatialPointsDataFrame(cart[,1:2], data.frame(cart[,3]))
  e <- extent(z)
  temp <- raster(z)
  res(temp) <- 10
  extent(temp) <- extent(-pi/2,pi/2,-pi/2,pi/2)
  ir <- rasterize(p, temp, fun = mean)
  ir <- ir$cart...3.
  extent(ir) <- e
  ir <- disaggregate(ir, 10, method = "bilinear")

  ir
}
