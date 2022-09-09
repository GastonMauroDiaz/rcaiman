#' Mask blue sky
#'
#' It is a wrapper function around \code{\link{membership_to_color}}. It
#' masks pixels that are likely blue sky.
#'
#' @inheritParams enhance_caim
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#'
#' @export
#'
#' @family Segmentation Functions
#'
#' @examples
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' m <- !is.na(z)
#' blue_sky <- mask_blue_sky(caim, m)
#' plot(blue_sky)
mask_blue_sky <- function(caim, m = NULL) {
  if (is.null(m)) m <- !is.na(caim$Red)
  mem <- membership_to_color(caim, sRGB(0.1,0.4,0.8))
  blue_sky <- mem$membership_to_target_color > 0.75
  blue_sky
}
