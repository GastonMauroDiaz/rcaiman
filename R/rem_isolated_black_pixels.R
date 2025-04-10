#' Remove isolated black pixels
#'
#' @param bin [SpatRaster-class]. Binarized canopy image.
#'
#' @returns An object of class [SpatRaster-class] with values `0` and `1`.
#' @export
#' @family Binarization Functions
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' sky <- read_ootb_sky_model(gsub(".txt", "", path), z, a)
#'
#' bin <- apply_thr(r/sky$sky, 0.9)
#' plot(bin)
#' bin2 <- rem_isolated_black_pixels(bin)
#' plot(bin2)
#' plot(bin2 - bin)
#' }
rem_isolated_black_pixels <- function(bin) {
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  ma <- matrix(c(1,1,1,1,-8,1,1,1,1), ncol = 3, nrow = 3)
  bin[terra::focal(bin, ma) == 8] <- 1
  as.logical(bin)
}
