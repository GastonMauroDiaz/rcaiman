#' Interpolate digital numbers
#'
#' Interpolate values from hemispherical photographs.
#'
#' This function use \code{\link[lidR]{knnidw}} as workhorse function, so
#' arguments \code{k}, \code{p}, and \code{rmax} are passed to it.
#'
#' This function is based on \insertCite{Lang2010;textual}{rcaiman}. In theory,
#' interpolation requires a linear relation between DNs and the amount of light
#' reaching the sensor. To that end, photographs should be taken in RAW format
#' to avoid gamma correction \insertCite{Lang2010}{rcaiman}. As a compromise
#' solution, \code{\link{gbc}} can be used.
#'
#' The vignetting effect also hindered linear relation between DNs and the
#' amount of light reaching the sensor. Please refer to
#' \insertCite{Lang2010;textual}{rcaiman} for more details about the vignetting
#' effect.
#'
#' The use of \code{k = 1} solve the linear dilemma from the theoretical point
#' of view since no averaging is taking place in the calculations.
#'
#' Default parameters are the used by \insertCite{Lang2010;textual}{rcaiman}.
#' The argument \code{rmax} should account for between 15 to 20 degrees, but is
#' expressed in pixels units. So, image resolution and lens projections should
#' be taken into account to properly set this argument.
#'
#'
#' @inheritParams fit_cie_sky_model
#' @inheritParams lidR::knnidw
#'
#' @references \insertRef{Lang2010}{rcaiman}
#'
#' @export
#'
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
#' bin <- ootb_mblt(blue, z, a)$bin
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_marks <- extract_sky_marks(blue, bin, g)
#' sky <- interpolate_dns(blue, sky_marks)
#' plot(sky)
#' }
interpolate_dns <- function(r, sky_marks,
                            k = 3,
                            p = 2,
                            rmax = 200,
                            use_window = TRUE) {
  .check_if_r_was_normalized(r)

  stopifnot(length(k) == 1)
  stopifnot(length(p) == 1)
  stopifnot(length(rmax) == 1)
  stopifnot(length(use_window) == 1)
  stopifnot(is.logical(use_window))
  stopifnot(is.numeric(k))
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(rmax))
  stopifnot(is.data.frame(sky_marks))

  cells <- cellFromRowCol(r, sky_marks$row, sky_marks$col)
  if (use_window) {
    xy <-  xyFromCell(r, cells)
    sky_marks$dn <-  extract(r, xy, buffer = 1.5, fun = "mean")
  } else {
    sky_marks$dn <-  r[cells]
  }

  las <- .make_fake_las(
    c(xy[, 1]     , 0 - rmax, 0 - rmax       , xmax(r) + rmax, xmax(r) + rmax),
    c(xy[, 2]     , 0 - rmax, ymax(r) + rmax , ymax(r) + rmax, 0 - rmax),
    c(sky_marks$dn, 0       , 0              ,              0, 0)
  )
  las@data$Classification <- 2
  ir <- lidR::grid_terrain(las, res = 1,
                           full_raster = TRUE,
                           algorithm = lidR::knnidw(k = k,
                                                    p = p,
                                                    rmax = rmax)
                           )
  resample(ir, r)
}



#' Interpolate by reprojecting
#'
#' Interpolate values from hemispherical photographs by reprojecting.
#'
#' This function assume an hemispherical image as input. It is a version of
#' \code{\link{interpolate_dns}} in which the interpolation takes place after
#' a cylindrical reprojection. To ensure hemispherical continuity, portions of
#' the reprojected image (30 degrees wide) are taken from the margin,
#' duplicated, translated, and mirrored as the task requires. The reprojected
#' image has resolution of one degree. Thus, \code{rmax} is in degrees.
#'
#'
#' @inheritParams ootb_mblt
#' @inheritParams zenith_image
#' @inheritParams fit_cie_sky_model
#' @inheritParams interpolate_dns
#'
#' @references \insertRef{Lang2010}{rcaiman}
#'
#' @noRd
#'
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
#' bin <- ootb_mblt(blue, z, a)$bin
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_marks <- extract_sky_marks(blue, bin, g)
#' sky <- interpolate_reproj(blue, z, a, lens("Nikon_FCE9"), sky_marks)
#' plot(sky)
#' }
interpolate_reproj <- function(r, z, a, lens_coef, sky_marks,
                            k = 3,
                            p = 2,
                            rmax = 20,
                            use_window = TRUE) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  stopifnot(ncol(r) == nrow(r))

  stopifnot(length(k) == 1)
  stopifnot(length(p) == 1)
  stopifnot(length(rmax) == 1)
  stopifnot(length(use_window) == 1)
  stopifnot(is.logical(use_window))
  stopifnot(is.numeric(k))
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(rmax))
  stopifnot(is.data.frame(sky_marks))

  stopifnot(is.numeric(lens_coef))

  cells <- cellFromRowCol(a, sky_marks$row, sky_marks$col)
  sky_marks$a <- a[cells]
  sky_marks$z <- z[cells]
  if (use_window) {
    xy <-  xyFromCell(r, cells)
    sky_marks$dn <-  extract(r, xy, buffer = 1.5, fun = "mean")
  } else {
    sky_marks$dn <-  r[cells]
  }

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
                                           full_raster = TRUE,
                                           algorithm = lidR::knnidw(k = k,
                                                                    p = p,
                                                                    rmax = rmax)
                                            )
                   ))
  ir <- crop(ir, extent(0,360,0,90))
  az <- xyFromCell(ir, 1:ncell(ir))
  rr <- calc_relative_radius(az[,2], lens_coef)
  pol <- data.frame(theta = az[,1] * pi/180 + pi/2,
                    r = rr * 90 * pi/180,
                    z = ir[])
  cart <- pracma::pol2cart(as.matrix(pol))
  p <- sp::SpatialPointsDataFrame(cart[,1:2], data.frame(cart[,3]))
  e <- extent(z)
  temp <- raster(z)
  res(temp) <- 10
  extent(temp) <- extent(-pi/2,pi/2,-pi/2,pi/2)
  ir <- rasterize(p, temp, fun = mean)
  ir <- ir$cart...3.
  extent(ir) <- e
  disaggregate(ir, 10, method = "bilinear")
}





