#' Interpolate sky points
#'
#' Interpolate values from canopy photographs.
#'
#' This function use [lidR::knnidw()] as workhorse function, so arguments `k`,
#' `p`, and `rmax` are passed to it.
#'
#' This function is based on \insertCite{Lang2010;textual}{rcaiman}. In theory,
#' the best result would be obtained with data showing a linear relation between
#' digital numbers and the amount of light reaching the sensor. See
#' [extract_radiometry()] and [read_caim_raw()] for further details. As a
#' compromise solution, [gbc()] can be used.
#'
#' Default parameters are the ones used by
#' \insertCite{Lang2010;textual}{rcaiman}. According to these authors, the
#' argument `rmax` should account for between 15 to 20 degrees, but it is
#' expressed in pixels units. So, image resolution and lens projections should
#' be taken into account to set this argument properly.
#'
#' @param sky_points An object of class *data.frame*. The data.frame returned by
#'   [extract_rel_radiance()] or [extract_dn()], or a
#'   *data.frame* with same basic structure and names.
#' @inheritParams extract_dn
#' @param k Numeric vector of length one. Number of k-nearest neighbors.
#' @param p Numeric vector of length one. Power for inverse-distance weighting.
#' @param rmax Numeric vector of length one. Maximum radius for searching
#'   k-nearest neighbors (knn).
#' @param col_id Numeric or character vector of length one. The number or name
#'   of the colum with the values to interpolate.
#'
#' @references \insertAllCited{}
#'
#' @return An object of class [SpatRaster-class].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' com <- compute_complementary_gradients(caim)
#' chroma <- max(com$blue_yellow, com$cyan_red)
#' bin <- apply_thr(chroma, thr_isodata(chroma[!is.na(chroma)]))
#' bin <- bin & apply_thr(com$blue_yellow, -0.2)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_black = 3)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' sky_points <- extract_dn(r, sky_points)
#'
#' sky <- interpolate_planar(sky_points, r, col_id = 3)
#' plot(sky)
#' plot(r/sky)
#' }
interpolate_planar <- function(sky_points,
                               r,
                               k = 3,
                               p = 2,
                               rmax = 200,
                               col_id = "rr") {
  .is_single_layer_raster(r)
  stopifnot(length(k) == 1)
  stopifnot(length(p) == 1)
  stopifnot(length(rmax) == 1)
  stopifnot(length(col_id) == 1)
  stopifnot(.is_whole(k))
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(rmax))
  stopifnot(is.data.frame(sky_points))

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  xy <- terra::xyFromCell(r, cells)

  if (max(sky_points[,col_id]) < 250) {
    const <- 10000
  } else {
    const <- 1
  }


  las <- .make_fake_las(
    c(xy[, 1]     , 0 - rmax, 0 - rmax       , xmax(r) + rmax, xmax(r) + rmax),
    c(xy[, 2]     , 0 - rmax, ymax(r) + rmax , ymax(r) + rmax, 0 - rmax),
    c(sky_points[,col_id] * const, 0       , 0              ,              0, 0)
  )
  las@data$Classification <- 2
  lidR::crs(las) <- 7589

  ir <- suppressWarnings(
        lidR::rasterize_terrain(las, res = 1,
                                algorithm = lidR::knnidw(k = k,
                                                         p = p,
                                                         rmax = rmax)
                                )
         )
  ir <- terra::resample(ir, r) / const


  p <- terra::vect(xy, "points")
  terra::crs(p) <- crs(r)

  buff <- terra::buffer(p, rmax)
  buff <- terra::rasterize(buff, r)
  ir[is.na(buff)] <- NA
  ir
}
