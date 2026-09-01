#' Extract digital numbers from image with sampling points
#'
#' @description
#' Obtain digital numbers from a raster at positions defined by points, with
#' optional local averaging.
#'
#' @details
#' Wraps [terra::extract()].
#'
#' @param r [terra::SpatRaster-class].
#' @param sampling_points `data.frame` with columns `"row"` and `"col"` (see
#'   [xy_from_row_col()]). Additional columns are ignored.
#' @param use_window logical vector of length one. If `TRUE` (default), the
#'   digital number at each point is the average of the values extracted
#'   from input `r` with a window of \eqn{3 \times 3} pixels centered on the
#'   pixel that includes the point. If `FALSE`, only the value of the central
#'   pixel is retrieved.
#'
#' @return `data.frame` containing the original `sampling_points` plus one column per
#'   layer in `r` (named after the layers).
#'
#' @note For instructions on manually digitizing points, see the “Digitizing
#' sky points with ImageJ” and “Digitizing sky points with QGIS” sections in
#' [fit_cie_model()].
#'
#' @seealso
#' [sample_sky_points()] for automatically sampling sky points.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(100, equidistant_coef)
#' a <- azimuth_image(z)
#' pts <- fibonacci_points(z, a, 30)
#' display_caim(is.na(z), sampling_points = pts)
#' plot(is.na(z))
#' points(pts$col, nrow(z) - pts$row, col = 3, pch = 4)
#' }
extract_dn <- function(r, sampling_points, use_window = TRUE) {
  .assert_spatraster(r)
  .check_columns(sampling_points, c("col", "row"), names(r))
  .check_vector(use_window, "logical", 1)

  if (nrow(sampling_points) > 1000 && use_window) {
    message(
      sprintf("The call with `use_window = TRUE` over %s points will be time consuming.",
              nrow(sampling_points))
    )
  }

  i <- cell_indices_from_col_row(sampling_points$row, sampling_points$col)

  cells <- terra::cellFromRowCol(r, i$row_index, i$col_index)
  xy <-  terra::xyFromCell(r, cells)
  if (use_window) {
    dn <- Map(function(x, y) {
      ma <- expand.grid(c(-1,0,1) + x, c(-1,0,1) + y)
      if (terra::nlyr(r) > 1) {
        return(apply(terra::extract(r, ma, method = "simple")[,-1], 2,
                     function(x) mean(x, na.rm = TRUE)))
      } else {
        return(mean(terra::extract(r, ma, method = "simple",
                                   na.rm = TRUE)[,-1]))
      }
    }, xy[,1], xy[,2]) %>% as.data.frame() %>% t %>% unname()
  } else {
    dn <- terra::extract(r, xy, method = "simple")[,]
  }
  sp_names <- c(colnames(sampling_points), names(r))
  sampling_points <- cbind(sampling_points, dn)
  colnames(sampling_points) <- sp_names
  sampling_points
}
