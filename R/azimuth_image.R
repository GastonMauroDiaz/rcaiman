#' Azimuth image
#'
#' Build a single-layer image with azimuth angles as pixel values, assuming
#' upwards-looking hemispherical photography with the optical axis vertically
#' aligned.
#'
#' @inheritParams ootb_mblt
#' @param orientation Azimuth angle (degrees) at which the top of the image was
#'   pointing at the moment of taking the picture. The usual field protocol is
#'   recording the angle at which the top of the camera points.
#'
#' @return An object of class [SpatRaster-class] with azimuth angles in degrees.
#'   If the `orientation` argument is zero, North (0º) is pointing up as in
#'   maps, but East (90º) and West (270º) are flipped regarding to maps. To
#'   understand why is that, do the following: take two flash-card size pieces
#'   of paper; put one on a table in front of you and draw on it a compass rose;
#'   take the other and hold it with your arms extended over your head and,
#'   following the directions of the compass rose in front of you, draw another
#'   one in the paper side that face down--It will be an awkward position, like
#'   if you were taking an upward-looking photo with a mobile device while
#'   looking at the screen--; finally, put it down and compare both compass
#'   roses.
#' @export
#'
#' @family Lens Functions
#'
#' @examples
#' z <- zenith_image(1490, lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' plot(a)
#' \dontrun{
#' a <- azimuth_image(z, 45)
#' plot(a)
#' }
azimuth_image <- function (z, orientation = 0)
{
  .is_single_layer_raster(z, "z")
  mask <- is.na(z)

  xy <- terra::xyFromCell(z, seq(length = ncell(z)))
  sph <- pracma::cart2sph(
    matrix(c(xy[, 1] - ncol(z) / 2,
             xy[, 2] - ncol(z) / 2,
             terra::values(z)), ncol = 3)
  )

  values(z) <- sph[, 1] * 180 / pi
  values(z) <- terra::values(abs(terra::trans(z) - 180))
  # above line is to orient North up and West left

  if (orientation != 0) {
    if (!requireNamespace("imager", quietly = TRUE)) {
      stop(paste("Package \"imager\" needed for this function to work.",
                 "Please install it."),
           call. = FALSE)
    }
    picture_cw_rotation <- orientation
    v <- imager::as.cimg(as.array(z)) %>% suppressWarnings()
    v <- imager::rotate_xy(v, -picture_cw_rotation,
                           ncol(z) / 2, nrow(z) / 2, interpolation = 0)
    terra::values(z) <- as.matrix(v)
  }

  z[mask] <- NA
  names(z) <- "Azimuth image"
  z
}
