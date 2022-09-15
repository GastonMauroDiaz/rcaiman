#' Interpolate sky points
#'
#' Interpolate values from canopy photographs.
#'
#' This function use \code{\link[lidR]{knnidw}} as workhorse function, so
#' arguments \code{k}, \code{p}, and \code{rmax} are passed to it.
#'
#' This method is based on \insertCite{Lang2010;textual}{rcaiman}. In theory,
#' interpolation requires a linear relation between DNs and the amount of light
#' reaching the sensor. To that end, photographs should be taken in RAW format
#' to avoid gamma correction \insertCite{Lang2010}{rcaiman}. As a compromise
#' solution, \code{\link{gbc}} can be used.
#'
#' The vignetting effect also hinders the linear relation between DNs and the
#' amount of light reaching the sensor. Please refer to
#' \insertCite{Lang2010;textual}{rcaiman} for more details about the vignetting
#' effect.
#'
#' The use of \code{k = 1} solves the linear dilemma from the theoretical point
#' of view since no averaging is taking place in the calculations. However,
#' probably, it is best to use \code{k} greater than 1.
#'
#' Default parameters are the ones used by
#' \insertCite{Lang2010;textual}{rcaiman}. The argument \code{rmax} should
#' account for between 15 to 20 degrees, but it is expressed in pixels units.
#' So, image resolution and lens projections should be taken into account to set
#' this argument properly.
#'
#' The argument \code{g} should be the same used to obtain \code{sky_points}.
#' The result will be limited to the cells with at least one pixel covered by
#' the convex hull of the sky points.
#'
#' @inheritParams extract_sky_points
#' @param sky_points An object of class \emph{data.frame}. The result of a call
#'   to \code{\link{extract_rl}} or \code{\link{extract_dn}}, or a
#'   \emph{data.frame} with same basic structure and names.
#' @param k Numeric vector of length one. Number of k-nearest neighbors.
#' @param p Numeric vector of length one. Power for inverse-distance weighting.
#' @param rmax Numeric vector of length one. Maximum radius where to search for
#'   knn.
#' @param col_id Numeric vector of length one. ID of the column with the values
#'   to interpolate.
#'
#' @references \insertAllCited{}
#'
#' @family Sky Reconstruction Functions
#'
#' @return An object of class \linkS4class{SpatRaster}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- ootb_obia(caim, z, a)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' r <- gbc(caim$Blue*255)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_rl(r, z, a, sky_points, NULL)
#' sky <- interpolate_sky_points(sky_points$sky_points, g)
#' plot(sky)
#'
#' #modify g if the goal is to get the whole sky
#' g <- !is.na(z)
#' sky <- interpolate_sky_points(sky_points$sky_points, g)
#' plot(sky)
#' plot(r/sky)
#'
#' #restricted view canopy photo
#' path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' plot(caim)
#' r <- gbc(caim$Blue)
#' caim <- normalize(caim)
#'
#' bin <- ootb_obia(caim)
#'
#' g <- chessboard(caim, 100)
#' plot(g)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_dn(r, sky_points)
#' head(sky_points)
#' sky <- interpolate_sky_points(sky_points, !is.na(r), col_id = 3)
#' plot(sky)
#' plot(r/sky)
#' }
interpolate_sky_points <- function(sky_points, g,
                                   k = 3,
                                   p = 2,
                                   rmax = 200,
                                   col_id = "rl") {
  .is_single_layer_raster(g)
  stopifnot(length(k) == 1)
  stopifnot(length(p) == 1)
  stopifnot(length(rmax) == 1)
  stopifnot(length(col_id) == 1)
  stopifnot(.is_integerish(k))
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(rmax))
  stopifnot(is.data.frame(sky_points))

  cells <- terra::cellFromRowCol(g, sky_points$row, sky_points$col)
  xy <- terra::xyFromCell(g, cells)
  las <- .make_fake_las(
    c(xy[, 1]     , 0 - rmax, 0 - rmax       , xmax(g) + rmax, xmax(g) + rmax),
    c(xy[, 2]     , 0 - rmax, ymax(g) + rmax , ymax(g) + rmax, 0 - rmax),
    c(sky_points[,col_id] * 10000, 0       , 0              ,              0, 0)
  )
  las@data$Classification <- 2
  # lidR::crs(las) <- sf::st_crs(g)
  lidR::crs(las) <- 7589

  ir <- suppressWarnings(
        lidR::rasterize_terrain(las, res = 1,
                                # shape =  "bbox",
                                algorithm = lidR::knnidw(k = k,
                                                         p = p,
                                                         rmax = rmax)
                                )
         )
  ir <- terra::resample(ir, g) / 10000

  indexes <- grDevices::chull(xy)
  p <- terra::vect(xy[indexes, ], "polygons")
  ds <- terra::extract(g, p)
  g <- terra::subst(g, unique(ds[,2]), 1)
  g <- g == 1

  ir[!g] <- NA
  ir
}
