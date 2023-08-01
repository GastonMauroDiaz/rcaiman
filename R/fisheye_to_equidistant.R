#' Fisheye to equidistant
#'
#' Fisheye to equidistant projection (also known as polar projection).
#'
#' The pixel values and their image coordinates are treated as points to be
#' reprojected and interpolated. To that end, this function use
#' \code{\link[lidR]{knnidw}} as workhorse function, so arguments \code{k},
#' \code{p}, and \code{rmax} are passed to it. If the user does not input values
#' to these arguments, both \code{k} and \code{p} are automatically defined by
#' default as follow. When a binarized image is provided as argument \code{r},
#' both parameters are set to \code{1}. Otherwise, they are set to \code{9} and
#' \code{2}, respectively.
#'
#' @param r \linkS4class{SpatRaster}.
#' @inheritParams ootb_mblt
#' @inheritParams interpolate_sky_points
#' @param rmax Numeric vector of length one. Maximum radius where to search for
#'   \emph{knn}. Increase this value if pixels with value \code{0} or
#'   \code{FALSE} appears where other values are expected.
#' @param radius Numeric integer of length one. Radius of the reprojected
#'   hemispherical image (i.e., the output). Default \code{NULL} is equivalent
#'   to input the radius of \code{r}.
#'
#' @export
#'
#' @family Lens Functions
#'
#' @examples
#' \donttest{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- apply_thr(caim$Blue, 0.5)
#' bin_equi <- fisheye_to_equidistant(bin, z, a)
#' plot(bin)
#' plot(bin_equi)
#' #use write_bin(bin, "path/file_name") to have a file ready
#' #for calculating LAI with CIMES, GLA, CAN-EYE, etc.
#'
#' #It can be used to reproject RGB photographs
#' plotRGB(caim*255)
#' caim <- fisheye_to_equidistant(caim, z, a)
#' plotRGB(caim*255)
#' }
fisheye_to_equidistant <- function(r, z, a,
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

    m <- !is.na(z) & !is.na(r)
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
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  r

}
