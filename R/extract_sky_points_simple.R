#' Extract sky points
#'
#' Extract sky points for model fitting
#'
#' This function will automatically sample sky pixels following this simple
#' strategy:
#' * mask the region of `r` above 15 and below 75 degrees of zenith angle,
#' * dividing the hemisphere into sectors of 15 degrees each
#' (see [sectors_segmentation()]),
#' * search for the maximum digital value in each sector (n = 24),
#' * dividing the hemisphere into rings of 5 degrees each
#' (see [rings_segmentation()],
#' * search for the maximum digital value in each ring (n = 12)
#' * combine these local maxima (n = 36).
#'
#' @inheritParams fit_trend_surface
#'
#' @return An object of the class *data.frame* with two columns named
#'   *col* and *row*.
#' @export
#'
#' @family Tool Functions
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' sky_points <- extract_sky_points_simple(r, z, a)
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
extract_sky_points_simple <- function(r, z, a) {
  m <- mask_hs(z, 15, 75)
  sky_points <- extract_sky_points(r, m,
                                   sectors_segmentation(a, 15))
  rbind(sky_points,
        extract_sky_points(r, m,
                           rings_segmentation(z, 5))
                           )
}
