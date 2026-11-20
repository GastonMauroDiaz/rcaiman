#' Remove sampling indices not selected by a mask
#'
#' @description
#' Retain sampling indices whose cells are selected by a logical mask.
#'
#' @details
#' Only cells where `m` is `TRUE` are retained. Cells marked `FALSE` or `NA`,
#' and indices outside the image, are excluded. Selection uses the indexed
#' cell without averaging neighboring values. Other columns do not affect
#' selection, including columns with missing values.
#'
#' Construct the mask to express the selection criterion. For example,
#' `!is.na(d$zenith)` selects cells with a zenith angle, including zenith and
#' nadir where azimuth is undefined. Use
#' `!is.na(g$off_axis) & !is.na(g$internal_azimuth)` to require both internal
#' angles. These conditions can be combined with a mask of eligible pixels.
#'
#' @param sampling_indices data frame with positive integer-like columns
#'   `col` and `row`, numbered from 1. Additional columns are preserved.
#' @param m logical [terra::SpatRaster-class] with one layer, on the grid
#'   to which `sampling_indices` refers. `TRUE` selects a cell; `FALSE` and
#'   `NA` exclude it.
#'
#' @return data frame containing the retained rows of `sampling_indices`,
#'   with the same columns and in their original order. No raster values are
#'   appended. Repeated indices are preserved.
#'
#' @seealso [add_directions()], [extract_dn()], [rem_nearby_indices()]
#' @export
#' @examples
#' r <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 3))
#' r[] <- c(1, NA, 3, 4, 5, 6)
#' sampling_indices <- data.frame(
#'   col = c(1, 2, 3), row = c(1, 1, 2), note = c(NA, "a", "b")
#' )
#' m <- !is.na(r)
#' rem_unselected(sampling_indices, m)
#' # Select by image values as well as availability.
#' rem_unselected(sampling_indices, m & r < 5)
rem_unselected <- function(sampling_indices, m) {
  .assert_single_layer(m)
  if (!terra::hasValues(m)) {
    stop("`m` must contain logical mask values.", call. = FALSE)
  }
  mask_values <- m[]
  # terra may represent an entirely missing logical mask as numeric NA.
  if (!is.logical(mask_values) && !all(is.na(mask_values))) {
    stop("`m` must contain logical mask values.", call. = FALSE)
  }
  .check_columns(sampling_indices, required_columns = c("col", "row"))
  .check_vector(sampling_indices$col, "integerish", sign = "positive",
                name = "sampling_indices$col")
  .check_vector(sampling_indices$row, "integerish", sign = "positive",
                name = "sampling_indices$row")
  if (nrow(sampling_indices) == 0) return(sampling_indices)
  cells <- terra::cellFromRowCol(m, sampling_indices$row, sampling_indices$col)
  selected <- terra::extract(m, cells)[[1]]
  sampling_indices[selected %in% TRUE, , drop = FALSE]
}
