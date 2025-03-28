#' Dysplay a canopy image
#'
#' Wrapper function for [EBImage::display()]
#'
#' @param caim [SpatRaster-class], such as the output of [read_caim()] or
#'   [read_caim_raw()], or a single layer raster such as the output of
#'   [apply_thr()].
#'
#' @return No return value. Called for side effects
#' @export
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' r <- normalize(caim$Blue)
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' display_caim(caim)
#' display_caim(bin)
#' display_caim(c(normalize(caim$Blue), bin))
#' }
display_caim <- function(caim) {
  .this_requires_EBImage()
  caim <- normalize(caim)
  as.array(caim) %>% EBImage::display()
}
