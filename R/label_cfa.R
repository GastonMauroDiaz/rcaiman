#' Label color filter array (CFA)
#'
#' @param r [SpatRaster-class].
#' @inheritParams read_caim
#'
#' @return [SpatRaster-class].
#' @noRd
label_cfa <- function(r,
                      upper_left = NULL,
                      width = NULL,
                      height = NULL) {
  .is_single_layer_raster(r)
  stopifnot(.is_even(ncol(r)/2))
  stopifnot(.is_even(nrow(r)/2))
  r[] <- c(rep(c(0,1), ncol(r)/2), rep(c(2,3), ncol(r)/2)) %>% rep(., nrow(r)/2)
  crop_caim(r, upper_left, width, height)
}
