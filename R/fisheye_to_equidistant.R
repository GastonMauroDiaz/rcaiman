#' Fisheye to equidistant
#'
#' Fisheye to equidistant projection (also known as polar projection).
#'
#' The pixel values and their image coordinates are treated as points to be
#' reprojected and interpolated. To that end, this function use [lidR::knnidw()]
#' as workhorse function, so arguments `k`, `p`, and `rmax` are passed to it.
#'
#' @param r [SpatRaster-class]. A fish-eye image.
#' @inheritParams calc_co
#' @inheritParams sky_grid_segmentation
#' @inheritParams interpolate_planar
#' @param rmax Numeric vector of length one. Maximum radius where to search for
#'   *knn*. Increase this value if pixels with value `0` or
#'   `FALSE` appears where other values are expected.
#' @param radius Numeric integer of length one. Radius of the reprojected
#'   hemispherical image (i.e., the output).
#'
#' @note Default value for the `radius` argument is equivalent to input the
#' radius of the `r` argument.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' calc_diameter(c(0.801, 0.178, -0.179), 1052/2, 86.2)
#' z <- zenith_image(1058,  c(0.801, 0.178, -0.179))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(532, 386)
#'
#' caim <- expand_noncircular(caim, z, zenith_colrow)
#' m <- !is.na(caim$Red) & select_sky_vault_region(z, 0, 86.2)
#' caim[!m] <- 0
#'
#' bin <- apply_thr(caim$Blue, thr_isodata(caim$Blue[m]))
#'
#' display_caim(caim$Blue, bin)
#'
#' caim <- gbc(caim, 2.2)
#' caim <- correct_vignetting(caim, z, c(-0.0546, -0.561, 0.22)) %>%
#'   normalize_minmax()
#'
#' caim2 <- fisheye_to_equidistant(caim$Blue, z, a, m, radius = 600)
#' bin2 <- fisheye_to_equidistant(bin, z, a, m, radius = 600)
#' bin2 <- apply_thr(bin2, 0.5) #to turn it logical
#' # Use write_bin(bin2, "path/file_name") to have a file ready
#' # to calcute LAI with CIMES, GLA, CAN-EYE, etc.
#'
#' m2 <- fisheye_to_equidistant(m, z, a, !is.na(z), radius = 600)
#'
#'
#' caim <- read_caim(path)
#' caim <- expand_noncircular(caim, z, zenith_colrow)
#' plotRGB(caim)
#' caim <- fisheye_to_equidistant(caim, z, a, !is.na(z), radius = 600)
#' caim[!m2] <- 0
#' plotRGB(caim)
#' }
fisheye_to_equidistant <- function(r, z, a, m,
                                   radius = NULL,
                                   k = 1,
                                   p = 1,
                                   rmax = 100)
  {
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(r, z)
  terra::compareGeom(r, a)

  if (is.null(radius)) radius <- ncol(r) / 2

  .fisheye_to_equidistant <- function(r) {
    new_r <- zenith_image(radius * 2, lens())
    terra::ext(new_r) <- terra::ext(-pi / 2, pi / 2, -pi / 2, pi / 2)

    terra::compareGeom(r, m)
    .is_single_layer_raster(m)
    .is_logic_and_NA_free(m)
    m <- !is.na(z) & !is.na(r) & m

    pol <- data.frame(theta = a[m] * pi / 180 + pi / 2,
                      r = z[m] * pi / 180,
                      z = r[m])
    names(pol) <- c("theta", "r", "z")
    cart <- pracma::pol2cart(as.matrix(pol))

    res <- terra::res(new_r)[1]

    const <- 10000

    res <- res * 1000
    las <- .make_fake_las(cart[,1]*1000, cart[,2]*1000, cart[,3]*const)
    las@data$Classification <- 2
    lidR::crs(las) <- 7589
    ir <- suppressWarnings(
      lidR::rasterize_terrain(las, res = res,
                              algorithm = lidR::knnidw(k = k,
                                                       p = p,
                                                       rmax = res * rmax)
      )
    )
    ir[is.na(ir)] <- 0
    ir / const
  }
  layer_names <- names(r)
  r <- Map(function(r) .fisheye_to_equidistant(r), as.list(r))
  r <- terra::rast(r)
  names(r) <- layer_names

  i <- terra::cellFromXY(r, matrix(c(0, 0), ncol = 2))
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  zenith_colrow <- terra::rowColFromCell(r, i) %>% as.numeric() %>% rev()

  z2 <- zenith_image(radius*2, lens())
  r <- expand_noncircular(r, z2, zenith_colrow)
  r[is.na(r)] <- 0
  r
}
