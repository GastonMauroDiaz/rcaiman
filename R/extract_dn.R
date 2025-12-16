#' Extract digital numbers from image with sampling points
#'
#' @description
#' Obtain digital numbers from a raster at positions defined by points, with
#' optional local averaging.
#'
#' @details Wraps [terra::extract()] to support a \eqn{3 \times 3} window
#' centered on each target pixel (local mean). When it is disabled, only the
#' central pixel value is retrieved.
#'
#' @param r [terra::SpatRaster-class]. Image from which the points were
#'   sampled (or any raster with identical dimensions).
#' @param sampling_points `data.frame` with columns `row` and `col` (raster
#'   coordinates).
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
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # See fit_cie_model() for details on below file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- extract_dn(caim, sky_points)
#' head(sky_points)
#'
#' # To aggregate DNs across points (excluding 'row' and 'col'):
#' apply(sky_points[, -(1:2)], 2, mean, na.rm = TRUE)
#' }
extract_dn <- function(r, sampling_points, use_window = TRUE) {
  .assert_spatraster(r)
  .check_sky_points(sampling_points)
  .check_vector(use_window, "logical", 1)

  cells <- terra::cellFromRowCol(r, sampling_points$row, sampling_points$col)
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
    # r_smooth <- terra::focal(r, 3, "mean")
    # dn <- terra::extract(r_smooth, xy, method = "simple")[,]
  } else {
    dn <- terra::extract(r, xy, method = "simple")[,]
  }
  sp_names <- c(colnames(sampling_points), names(r))
  sampling_points <- cbind(sampling_points, dn)
  colnames(sampling_points) <- sp_names
  sampling_points
}
