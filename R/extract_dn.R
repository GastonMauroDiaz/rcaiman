#' Extract digital numbers from cell indices
#'
#' @description
#' Obtain digital numbers from a raster at positions defined by points, with
#' optional local averaging.
#'
#' @details
#' Extract raster values using [terra::extract()]. When `use_window = TRUE`,
#' use each indexed cell as the center of the extraction window.
#'
#' @param r [terra::SpatRaster-class].
#' @param sampling_indices `data.frame` with positive integer-like columns
#'   `"col"` and `"row"`, representing cell indices numbered from 1.
#'   Additional columns are preserved.
#' @param use_window logical vector of length one. If `TRUE` (default), the
#'   digital number at each point is the average of the values extracted
#'   from input `r` with a window of \eqn{3 \times 3} pixels centered on the
#'   pixel that includes the point. If `FALSE`, only the value of the central
#'   pixel is retrieved. Return `NA` if the window includes cells outside `r`
#'   or cells with missing values.
#'
#' @note Rasters without values return `NA` for every sampling index and layer.
#'   An input with zero rows returns a zero-row `data.frame` with the expected
#'   output columns.
#'
#' @return `data.frame` containing the original `sampling_indices` plus one
#'   column per layer in `r`, named after the layers.
#'
#' @seealso
#' [compute_features()], [rem_nearby_indices()], [add_directions()].
#'
#' @export
#'
#' @examples
#' r <- terra::rast(nrows = 3, ncols = 3)
#' r <- set_rcaiman_geometry(r)
#' r[] <- 1:9
#' sampling_indices <- data.frame(
#'   col = c(1, 2),
#'   row = c(1, 2)
#' )
#' extract_dn(r, sampling_indices, use_window = FALSE)
extract_dn <- function(r, sampling_indices, use_window = TRUE) {
  .assert_spatraster(r)
  .check_columns(sampling_indices, c("col", "row"), names(r))
  .check_vector(
    sampling_indices$col,
    "integerish",
    sign = "positive",
    name = "sampling_indices$col"
  )
  .check_vector(
    sampling_indices$row,
    "integerish",
    sign = "positive",
    name = "sampling_indices$row"
  )
  .check_vector(use_window, "logical", 1)

  if (nrow(sampling_indices) == 0 || !terra::hasValues(r)) {
    dn <- matrix(
      NA_real_,
      nrow = nrow(sampling_indices),
      ncol = terra::nlyr(r)
    )
    colnames(dn) <- names(r)
    return(cbind(sampling_indices, dn))
  }

  if (nrow(sampling_indices) > 1000 && use_window) {
    message(
      sprintf(
        "The call with `use_window = TRUE` over %s points will be time consuming.",
        nrow(sampling_indices)
      )
    )
  }

  cells <- terra::cellFromRowCol(
    r,
    sampling_indices$row,
    sampling_indices$col
  )
  xy <-  terra::xyFromCell(r, cells)
  if (use_window) {
    dn <- Map(function(x, y) {
      ma <- expand.grid(c(-1,0,1) + x, c(-1,0,1) + y)
      if (terra::nlyr(r) > 1) {
        return(apply(terra::extract(r, ma, method = "simple")[,-1], 2,
                     mean))
      } else {
        return(mean(terra::extract(r, ma, method = "simple",
                                   na.rm = TRUE)[,-1]))
      }
    }, xy[,1], xy[,2]) %>% as.data.frame() %>% t %>% unname()
  } else {
    dn <- terra::extract(r, xy, method = "simple")[,]
  }
  sp_names <- c(colnames(sampling_indices), names(r))
  sampling_indices <- cbind(sampling_indices, dn)
  colnames(sampling_indices) <- sp_names
  sampling_indices
}
