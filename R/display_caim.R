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
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' r <- normalize_minmax(caim$Blue)
#' g <- rings_segmentation(z, 30)
#' bin <- regional_thresholding(r, g,
#'                              method = "thr_isodata")
#' display_caim(caim, bin, g)
#' }
display_caim <- function(caim = NULL, bin = NULL, g = NULL) {
  .this_requires_EBImage()
  if (!is.null(g)) {
    laplacian <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), nrow = 3)
    g <- terra::focal(g, laplacian)
    g <- g != 0
  }
  if (!is.null(caim)) {
    x <- normalize_minmax(caim)
    if (!is.null(bin) & !is.null(g)) {
      x <- c(x, bin, g)
    } else {
      if(!is.null(bin)) {
        x <- c(x, bin)
      }  else if (!is.null(g)) {
        x <- c(x, g)
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
    x <- terra::t(x) #needed for compatibility
    as.array(x) %>% EBImage::display()
  } else {
    warning("Nothing to display")
  }
}
