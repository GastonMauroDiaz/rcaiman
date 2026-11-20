#' Crop and reposition a canopy image
#'
#' @description
#' Extract a rectangular region of cells from an image or its angular geometry,
#' or reposition retained cells in their original grid.
#'
#' @param r numeric [terra::SpatRaster-class] with one or more layers.
#' @param upper_left_indices positive integer-like numeric vector of length two.
#'   Cell indices of the upper-left cell included in the crop, in the order
#'   `c(col, row)`. See [y_down_from_y()].
#' @param width,height positive integer-like numeric vectors of length one.
#'   Number of columns and rows included in the crop.
#' @param original_dim positive integer-like numeric vector of length two.
#'   Original grid dimensions in the order `c(width, height)`.
#'
#' @details
#' Crops are aligned with cell boundaries and never divide or resample a cell.
#' The requested rectangle must lie entirely within the input image. For
#' [uncrop_caim()], supply the same `upper_left_indices` and the original grid
#' dimensions. The retained cells must fit entirely within that grid.
#'
#' Apply the same selection to an image and its internal geometry to keep their
#' cells aligned. Angular values remain valid after cropping and repositioning.
#'
#' The result is assigned the standardized raster geometry described in
#' [set_rcaiman_geometry()]. Continuous image coordinates change with this new
#' frame; use [crop_points()] to translate marks accordingly. [uncrop_caim()]
#' places retained cells back in the original frame, filling the remaining cells
#' with `NA`; use [uncrop_points()] to translate marks accordingly. Values
#' discarded by cropping are not recovered.
#'
#' @return
#' - `crop_caim()` returns a [terra::SpatRaster-class] containing the selected
#'   cells, with the same layers and values as `r`.
#' - `uncrop_caim()` returns a [terra::SpatRaster-class] with the original grid
#'   dimensions and the same layers as `r`. Retained cells preserve their
#'   values; all other cells contain `NA`.
#'
#' @seealso [read_caim()], [set_internal_geometry()], [crop_points()],
#'   [uncrop_points()]
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(30, 20),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(60, 40), geometry)
#'
#' caim <- g$off_axis + 20
#' original_dim <- c(ncol(caim), nrow(caim))
#' upper_left_indices <- c(2, 3)
#' cropped_caim <- crop_caim(
#'   caim,
#'   upper_left_indices = upper_left_indices,
#'   width = 10,
#'   height = 10
#' )
#' plot(cropped_caim)
#' dim(cropped_caim)
#' restored_caim <- uncrop_caim(cropped_caim, upper_left_indices, original_dim)
#' dim(restored_caim)
#' plot(restored_caim) # only the frame is restored
#'
#' # Use the same frame for the internal geometry.
#' geometry <- list(
#'   center_of_vision = floor(original_dim / 2),
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   normalization_radius = min(original_dim) / 2,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(original_dim, geometry)
#' g_crop <- crop_caim(g, upper_left_indices, width = 5, height = 4)
#' g_back <- uncrop_caim(g_crop, upper_left_indices, original_dim)
#'
#' # Define a square around the center of vision.
#' center_of_vision <- c(5, 6)
#' crop_radius <- 2
#' # The center uses boundary coordinates; cell indices start at 1.
#' upper_left_indices <- center_of_vision - crop_radius + 1
#' side <- 2 * crop_radius
#' centered_caim <- crop_caim(caim, upper_left_indices, side, side)
#' plot(centered_caim)
crop_caim <- function(r, upper_left_indices, width, height) {
  .assert_spatraster(r)
  .check_vector(upper_left_indices, "integerish", 2, sign = "positive")
  .check_vector(width, "integerish", 1, sign = "positive")
  .check_vector(height, "integerish", 1, sign = "positive")

  first_col <- upper_left_indices[1]
  first_row <- upper_left_indices[2]
  last_col <- first_col + width - 1
  last_row <- first_row + height - 1

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

  r <- terra::crop(set_rcaiman_geometry(r), crop_extent)
  set_rcaiman_geometry(r)
}
