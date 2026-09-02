#' Crop a canopy image by cell indices
#'
#' @description
#' Extract an exact rectangular region of cells from a canopy image. The
#' `upper_left_indices` argument identifies the first cell included in the
#' result, and `width` and `height` give the number of included columns and
#' rows.
#'
#' @param r [terra::SpatRaster-class].
#' @param upper_left_indices integer-like numeric vector of length two. Cell
#'   indices of the upper-left cell included in the crop, in the order
#'   `c(col, row)`. Cell indices are numbered from 1. If `NULL`, `width` and
#'   `height` must also be `NULL`, and no crop is applied.
#' @param width,height positive integer-like numeric vectors of length one.
#'   Number of columns and rows included in the crop. If `NULL`,
#'   `upper_left_indices` must also be `NULL`, and no crop is applied.
#'
#' @details
#' Cropping is aligned with cell boundaries and never divides or resamples a
#' cell. For example, `upper_left_indices = c(4, 3)`, `width = 2`, and
#' `height = 3` retain columns 4--5 and rows 3--5.
#'
#' The result is assigned the standardized raster geometry described in
#' [set_rcaiman_geometry()].
#'
#' @return [terra::SpatRaster-class] containing the selected cells, with the
#'   same layers and values as `r`.
#'
#' @export
#'
#' @examples
#' caim <- read_caim()
#' cropped_caim <- crop_caim(
#'   caim,
#'   upper_left_indices = c(2, 3),
#'   width = 5,
#'   height = 4
#' )
#' dim(cropped_caim)
crop_caim <- function(r,
                      upper_left_indices = NULL,
                      width = NULL,
                      height = NULL) {
  .assert_spatraster(r)
  .check_vector(
    upper_left_indices,
    "integerish",
    2,
    allow_null = TRUE,
    sign = "positive"
  )
  .check_vector(
    width,
    "integerish",
    1,
    allow_null = TRUE,
    sign = "positive"
  )
  .check_vector(
    height,
    "integerish",
    1,
    allow_null = TRUE,
    sign = "positive"
  )

  supplied <- !vapply(
    list(upper_left_indices, width, height),
    is.null,
    logical(1)
  )
  if (any(supplied) && !all(supplied)) {
    stop(
      "`upper_left_indices`, `width`, and `height` must be supplied together.",
      call. = FALSE
    )
  }

  r <- set_rcaiman_geometry(r)
  if (!all(supplied)) {
    return(r)
  }

  first_col <- as.integer(upper_left_indices[1])
  first_row <- as.integer(upper_left_indices[2])
  last_col <- first_col + as.integer(width) - 1L
  last_row <- first_row + as.integer(height) - 1L

  if (last_col > terra::ncol(r) || last_row > terra::nrow(r)) {
    stop(
      "The crop defined by `upper_left_indices`, `width`, and `height` ",
      "extends outside `r`.",
      call. = FALSE
    )
  }

  x_min <- first_col - 1
  x_max <- last_col
  y_max <- terra::nrow(r) - first_row + 1
  y_min <- terra::nrow(r) - last_row
  crop_extent <- terra::ext(x_min, x_max, y_min, y_max)

  r <- terra::crop(r, crop_extent)
  set_rcaiman_geometry(r)
}
