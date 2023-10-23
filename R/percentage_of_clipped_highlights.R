#' Percentage of clipped highlights
#'
#' Wrapper function for [terra::freq()]
#'
#' @param r Single-layer object from the [SpatRaster-class].
#' @inheritParams masking
#'
#' @return Numeric vector of lenght one.
#' @export
#'
#' @family Tool Functions
#'
#' @examples
#' r <- read_caim()$Blue
#' z <- zenith_image(ncol(r), lens())
#' m <- !is.na(z)
#' percentage_of_clipped_highlights(r, m)
#' r <- normalize(r, 0, 1000, TRUE)
#' percentage_of_clipped_highlights(r, m)
percentage_of_clipped_highlights <- function(r, m) {
  .is_single_layer_raster(r)
  .was_normalized(r)
  .is_logic_and_NA_free(m)
  terra::compareGeom(r, m)
  (terra::freq(r == 1, value = 1)$count/sum(m[]) * 100) %>% round()
}
