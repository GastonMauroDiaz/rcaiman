.crop_offset_coord <- function(
    upper_left_indices = NULL,
    snapped_zenith_coord = NULL,
    horizon_radius = NULL,
    allow_null = FALSE
) {
  .check_vector(
    upper_left_indices,
    "integerish",
    2,
    allow_null = TRUE,
    sign = "positive"
  )
  .check_vector(
    snapped_zenith_coord,
    "integerish",
    2,
    allow_null = TRUE,
    sign = "positive"
  )
  .check_vector(
    horizon_radius,
    "integerish",
    1,
    allow_null = TRUE,
    sign = "positive"
  )

  free_supplied <- !is.null(upper_left_indices)
  centered_supplied <- c(
    !is.null(snapped_zenith_coord),
    !is.null(horizon_radius)
  )

  if (free_supplied && any(centered_supplied)) {
    stop(
      "Supply either `upper_left_indices` or the centered specification, ",
      "not both.",
      call. = FALSE
    )
  }
  if (any(centered_supplied) && !all(centered_supplied)) {
    stop(
      "`snapped_zenith_coord` and `horizon_radius` must be supplied together.",
      call. = FALSE
    )
  }
  if (!free_supplied && !any(centered_supplied)) {
    if (allow_null) {
      return(NULL)
    }
    stop(
      "Supply `upper_left_indices` or both `snapped_zenith_coord` and ",
      "`horizon_radius`.",
      call. = FALSE
    )
  }

  if (free_supplied) {
    return(upper_left_indices - 1)
  }
  snapped_zenith_coord - horizon_radius
}


#' Crop a canopy image
#'
#' @description
#' Extract an exact region of cells from a canopy image. Define either a free
#' rectangular crop or a square crop centered on a snapped zenith.
#'
#' @param r [terra::SpatRaster-class].
#' @param upper_left_indices integer-like numeric vector of length two. Cell
#'   indices of the upper-left cell included in the crop, in the order
#'   `c(col, row)`. Cell indices are numbered from 1. If `NULL`, `width` and
#'   `height` must also be `NULL`.
#' @param width,height positive integer-like numeric vectors of length one.
#'   Number of columns and rows included in the crop. If `NULL`,
#'   `upper_left_indices` must also be `NULL`.
#' @param snapped_zenith_coord integer-like numeric vector of length two.
#'   Continuous image coordinates of the optical center in the input image
#'   frame, in the order `c(x, y_down)`. The coordinates identify the center of
#'   the nearest 2 x 2 cell window. If `NULL`, `horizon_radius` must also be
#'   `NULL`.
#' @param horizon_radius positive integer-like numeric vector of length one.
#'   Radius of the declared processing geometry in cell units. If `NULL`,
#'   `snapped_zenith_coord` must also be `NULL`.
#'
#' @details
#' Crops are aligned with cell boundaries and never divide or resample a cell.
#' Supply either `upper_left_indices`, `width`, and `height`, or
#' `snapped_zenith_coord` and `horizon_radius`. The two specifications are
#' mutually exclusive. If all five arguments are `NULL`, no crop is applied.
#'
#' For a free crop, `upper_left_indices = c(4, 3)`, `width = 2`, and
#' `height = 3` retain columns 4--5 and rows 3--5. For a centered crop, the
#' first retained cell is identified by
#' `snapped_zenith_coord - horizon_radius + 1`, and both sides have length
#' `2 * horizon_radius`.
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
#'
#' centered_caim <- crop_caim(
#'   caim,
#'   snapped_zenith_coord = c(5, 6),
#'   horizon_radius = 2
#' )
crop_caim <- function(r,
                      upper_left_indices = NULL,
                      width = NULL,
                      height = NULL,
                      snapped_zenith_coord = NULL,
                      horizon_radius = NULL) {
  .assert_spatraster(r)
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

  rectangle_supplied <- !vapply(
    list(upper_left_indices, width, height),
    is.null,
    logical(1)
  )
  centered_supplied <- !vapply(
    list(snapped_zenith_coord, horizon_radius),
    is.null,
    logical(1)
  )

  if (any(rectangle_supplied) && any(centered_supplied)) {
    stop(
      "Supply either the free or centered crop specification, not both.",
      call. = FALSE
    )
  }
  if (any(rectangle_supplied) && !all(rectangle_supplied)) {
    stop(
      "`upper_left_indices`, `width`, and `height` must be supplied together.",
      call. = FALSE
    )
  }

  offset_coord <- .crop_offset_coord(
    upper_left_indices,
    snapped_zenith_coord,
    horizon_radius,
    allow_null = TRUE
  )

  r <- set_rcaiman_geometry(r)
  if (is.null(offset_coord)) {
    return(r)
  }

  if (all(centered_supplied)) {
    width <- height <- horizon_radius * 2
  }

  first_col <- as.integer(offset_coord[1] + 1)
  first_row <- as.integer(offset_coord[2] + 1)
  last_col <- first_col + as.integer(width) - 1L
  last_row <- first_row + as.integer(height) - 1L

  outside_r <- first_col < 1L || first_row < 1L ||
    last_col > terra::ncol(r) || last_row > terra::nrow(r)
  if (outside_r) {
    stop(
      "The requested crop extends outside `r`.",
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
