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
#' @param sky_points An object of class \emph{data.frame}. The result of a call
#'   to \code{\link{extract_rl}} or \code{\link{extract_dn}}, or a
#'   \emph{data.frame} with same basic structure and names.
#' @param r \linkS4class{SpatRaster}. The image from which \code{sky_points}
#'   was obtained.
#' @param k Numeric vector of length one. Number of k-nearest neighbors.
#' @param p Numeric vector of length one. Power for inverse-distance weighting.
#' @param rmax Numeric vector of length one. Maximum radius where to search for
#'   \emph{knn}.
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
#' \donttest{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' bin <- ootb_obia(caim, z, a)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' r <- gbc(caim$Blue*255)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_rl(r, z, a, sky_points, NULL)
#' sky <- interpolate_sky_points(sky_points$sky_points, r)
#' plot(sky)
#'
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
#' sky <- interpolate_sky_points(sky_points, r, col_id = 3)
#' plot(sky)
#' plot(r/sky)
#' }
interpolate_sky_points <- function(sky_points,
                                   r,
                                   k = 3,
                                   p = 2,
                                   rmax = 200,
                                   col_id = "rl") {
  .is_single_layer_raster(r)
  stopifnot(length(k) == 1)
  stopifnot(length(p) == 1)
  stopifnot(length(rmax) == 1)
  stopifnot(length(col_id) == 1)
  stopifnot(.is_integerish(k))
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(rmax))
  stopifnot(is.data.frame(sky_points))

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  xy <- terra::xyFromCell(r, cells)
  las <- .make_fake_las(
    c(xy[, 1]     , 0 - rmax, 0 - rmax       , xmax(r) + rmax, xmax(r) + rmax),
    c(xy[, 2]     , 0 - rmax, ymax(r) + rmax , ymax(r) + rmax, 0 - rmax),
    c(sky_points[,col_id] * 10000, 0       , 0              ,              0, 0)
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
  ir <- terra::resample(ir, r) / 10000


  p <- terra::vect(xy, "points")
  terra::crs(p) <- crs(r)

  buff <- terra::buffer(p, rmax)
  buff <- terra::rasterize(buff, r)
  ir[is.na(buff)] <- NA
  ir
}
