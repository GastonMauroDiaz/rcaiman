#' Build CFA index image
#'
#' @description
#' Index pixels according to their position in a 2 x 2 CFA block.
#'
#' @details
#' The raster dimensions must be even to accommodate an integer number of
#' 2 x 2 CFA blocks. Labels are assigned as:
#' \deqn{
#' \begin{matrix}
#' 0 & 1 \\
#' 2 & 3
#' \end{matrix}
#' }
#'
#' The input values are ignored. Only the raster dimensions are used.
#'
#'
#' @param r single-layer [terra::SpatRaster-class].
#'
#' @return Integer [terra::SpatRaster-class] with the same dimensions as `r`
#'   and values `0`, `1`, `2`, and `3`.
#'
#' @export
cfa_index_image <- function(r) {
  .assert_single_layer(r)
  .check_vector(ncol(r), "even_integerish", 1)
  .check_vector(nrow(r), "even_integerish", 1)
  r[] <- c(rep(c(0,1), ncol(r)/2), rep(c(2,3), ncol(r)/2)) %>% rep(., nrow(r)/2)
  r
}
