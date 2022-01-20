#' Azimuth image
#'
#' Build a single layer image with azimuth angles as pixel values.
#'
#' @param z \linkS4class{RasterLayer} built with
#'   \code{\link{zenith_image}}.
#'
#'
#' @return An object of class \linkS4class{RasterLayer} of azimuth angles
#'   in degrees. North (0º) is pointing up as in maps, but East (90º) and West
#'   (270º) are flipped respect to maps. To understand why is that, take two
#'   flash-card size pieces of paper. Put one on a table in front of you and
#'   draw on it a compass rose. Take the other and hold it with your arms
#'   extended over your head, and, following the directions of the compass rose
#'   in front of you, draw another compass rose in the paper side that face
#'   down. Then, put it down and compare both compass roses.
#' @export
#'
#' @family Lens functions
#'
#' @examples
#' z <- zenith_image(1490, lens("Nikon_FCE9"))
#' azimuth_image(z)
#' plot(z)
azimuth_image <- function (z)
{
  stopifnot(class(z) == "RasterLayer")
  mask <- is.na(z)

  xy <- xyFromCell(z, seq(length = ncell(z)))
  v <- values(z)
  sph <- pracma::cart2sph(
    matrix(c(xy[, 1] - ncol(z) / 2, xy[, 2] - ncol(z) / 2, values(z)), ncol = 3)
  )

  values(z) <- sph[, 1] * 180 / pi
  values(z) <- values(abs(raster::t(z) - 180)) #to orient North up and West left

  z[mask] <- NA

  z
}
