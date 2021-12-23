#' Azimuth image
#'
#' Build a single layer image with azimuth angles as pixel values.
#'
#' @param z \code{\linkS4class{RasterLayer}} built with
#'   \code{\link{zenith_image}}.
#'
#'
#' @return \code{\linkS4class{RasterLayer}}.
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
