#' Optimize distance to black
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams extract_sky_points
#' @param bin [SpatRaster-class]. This should be a preliminary binarization of
#'   `r` useful for masking pixels that are very likely pure sky pixels.
#' @inheritParams calc_co
#'
#'
#'
#' @returns Numeric vector of length one.
#' @export
#'
optim_dist_to_black <- function(r, z, a, m, bin, g) {
  g30 <- sky_grid_segmentation(z, a, 30)
  g30[!m] <- 0
  dist_to_black <- 11
  sampling_pct <- 0
  while (sampling_pct < 100 & dist_to_black > 3) {
    dist_to_black <- dist_to_black - 2
    sky_points <- extract_sky_points(r, bin, g,
                                     dist_to_black = dist_to_black)
    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
      (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
  }
  if (sampling_pct < 75) {
    dist_to_black <- 1
    sky_points <- extract_sky_points(r, bin, g,
                                     dist_to_black = dist_to_black)
    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
      (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
  }
  if (sampling_pct < 50) {
    dist_to_black <- NULL
  }
  dist_to_black
}
