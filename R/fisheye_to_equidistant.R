#' Fisheye to equidistant
#'
#' Fisheye to equidistant projection (also known as polar projection).
#'
#' The pixel values and their image coordinates are treated as points to be
#' reprojected and interpolated. To that end, this function use [lidR::knnidw()]
#' as workhorse function, so arguments `k`, `p`, and `rmax` are passed to it. If
#' the user does not input values to these arguments, both `k` and `p` are
#' automatically defined by default as follow: when a binarized image is
#' provided as argument `r`, both parameters are set to `1`; otherwise, they are
#' set to `9` and `2`, respectively.
#'
#' @param r [SpatRaster-class]. A fish-eye image.
#' @inheritParams calc_co
#' @inheritParams ootb_mblt
#' @inheritParams interpolate_sky_points
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
#' @family Lens Functions
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' r <- correct_vignetting(r, z, c(0.0638, -0.101)) %>% normalize()
#' bin <- ootb_mblt(r, z, a)$bin
#' bin_equi <- fisheye_to_equidistant(bin, z, a)
#' plot(bin)
#' plot(bin_equi)
#' # Use write_bin(bin, "path/file_name") to have a file ready
#' # to calcute LAI with CIMES, GLA, CAN-EYE, etc.
#'
#' # It can be used to reproject RGB photographs
#' plotRGB(caim)
#' caim <- fisheye_to_equidistant(caim, z, a)
#' plotRGB(caim)
#' }
fisheye_to_equidistant <- function(r, z, a,
                                   m = NULL,
                                   radius = NULL,
                                   k = NULL,
                                   p = NULL,
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

    if(!is.null(m)) {
      terra::compareGeom(r, m)
      .is_single_layer_raster(m)
      .is_logic_and_NA_free(m)
      m <- !is.na(z) & !is.na(r) & m
    } else {
      m <- !is.na(z) & !is.na(r)
    }
    pol <- data.frame(theta = a[m] * pi / 180 + pi / 2,
                      r = z[m] * pi / 180,
                      z = r[m])
    names(pol) <- c("theta", "r", "z")
    cart <- pracma::pol2cart(as.matrix(pol))

    res <- terra::res(new_r)[1]

    if (is.logical(r[1][,])) {
      if (is.null(k)) k <- 1
      if (is.null(p)) p <- 1

      las <- .make_fake_las(cart[,1], cart[,2], cart[,3])
      las@data$Classification <- 2
      lidR::crs(las) <- 7589

      ir <- suppressWarnings(
        lidR::rasterize_terrain(las, res = res,
                                algorithm = lidR::knnidw(k = k,
                                                         p = p,
                                                         rmax = res * rmax)
        )
      )
      terra::ext(ir) <- terra::ext(0, ncol(ir), 0, nrow(ir))
      ir <- apply_thr(ir, 0.5)
    } else {
      if (is.null(k)) k <- 9
      if (is.null(p)) p <- 2
      if (max(cart[,3]) < 1000) {
        const <- 10000
      } else {
        const <- 1
      }

      las <- .make_fake_las(cart[,1], cart[,2], cart[,3]*const)
      las@data$Classification <- 2
      lidR::crs(las) <- 7589

      ir <- suppressWarnings(
        lidR::rasterize_terrain(las, res = res,
                                algorithm = lidR::knnidw(k = k,
                                                         p = p,
                                                         rmax = res * rmax)
        )
      )
      terra::ext(ir) <- terra::ext(0, ncol(ir), 0, nrow(ir))
      ir[is.na(ir)] <- 0
      ir <- ir / const
    }
    ir
  }

  if (terra::nlyr(r) == 1) {
    r <- suppressWarnings(.fisheye_to_equidistant(r))
  } else {
    layer_names <- names(r)
    r <- Map(function(r) .fisheye_to_equidistant(r), as.list(r))
    r <- terra::rast(r)
    names(r) <- layer_names
  }
  terra::ext(r) <- terra::ext(0, radius*2, 0, radius*2)
  if(terra::res(r)[1] != 1) {
    template <- rast(r)
    terra::res(template) <- 1
    r <- terra::resample(r, template, method = "near")
  }
  r

}
