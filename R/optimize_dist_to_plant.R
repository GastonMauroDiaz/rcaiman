#' Optimize the 'dist_to_plant' argument
#'
#' Find the distance to the class plant that balance pure sky retrival with
#' sampling as much sky as possible.
#'
#' @inheritParams ootb_mblt
#' @param m [SpatRaster-class]. A mask, check [mask_hs()].
#' @param bin
#'
#' @return A numeric vector of length one.
#' @family Tool Functions
#' @seealso [ootb_sky_reconstruction()]
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r <- normalize(caim$Blue)
#'
#' bin <- find_sky_pixels(r, z, a)
#' bin <- ootb_mblt(r, z, a, bin)
#' plot(bin$bin)
#'
#' dist_to_plant <- optimize_dist_to_plan(r, z, a, m, bin)
#' ## use this parameter on ootb_sky_reconstruction()
#' }
optimize_dist_to_plant <- function(r, z, a, m, bin) {
  g30 <- sky_grid_segmentation(z, a, 30)
  g30[!m] <- 0
  g <- sky_grid_segmentation(z, a, 10)
  dist_to_plant <- 11
  sampling_pct <- 0
  while (sampling_pct < 100 & dist_to_plant > 3) {
    dist_to_plant <- dist_to_plant-2
    sky_points <- extract_sky_points(r, bin, g, dist_to_plant = dist_to_plant)
    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
      (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
  }
  if (sampling_pct < 75) {
    dist_to_plant <- 1
  }
  dist_to_plant
}
