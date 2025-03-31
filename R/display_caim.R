#' Dysplay a canopy image
#'
#' Wrapper function for [EBImage::display()]
#'
#' @param caim [SpatRaster-class], such as the output of [read_caim()] or
#'   [read_caim_raw()], or a single layer raster such as the output of
#'   [enhance_caim()].
#' @param bin description
#' @param g description
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
display_caim <- function(caim = NULL, bin = NULL, g = NULL) {
  .this_requires_EBImage()
  if (!is.null(g)) {
    g <- terra::focal(g, 3, sd)
    g <- g != 0
  }
  if (!is.null(caim)) {
    x <- normalize(caim)
    if (!is.null(bin) & !is.null(g)) {
      x <- c(x, bin, g)
    } else {
      if(!is.null(bin)) {
        x <- c(x, bin)
      }  else if (!is.null(g)) {
        x <- c(x, bin)
      }
    }
  } else {
    if (!is.null(bin) & !is.null(g)) {
      x <- c(bin, g)
    } else {
      if(!is.null(bin)) {
        x <- bin
      }  else {
        x <- g
      }
    }
  }
  if (!is.null(x)) {
    as.array(x) %>% EBImage::display()
  } else {
    warning("Nothing to display")
  }
}
