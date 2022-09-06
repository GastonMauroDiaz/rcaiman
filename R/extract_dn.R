#' Extract digital numbers
#'
#' It is a wrapper function around \code{\link[terra]{extract}}.
#'
#' @inheritParams fisheye_to_equidistant
#' @inheritParams extract_rl
#'
#' @return An object of the class \emph{data.frame} which is the argument
#'   \code{sky_points} with an added column for each layer from \code{r}. Layer names are
#'   used to name the new columns.
#' @export
#'
#' @examples
#' caim <- read_caim()
#' r <- gbc(caim$Blue)
#' bin <- apply_thr(r, thr_isodata(r[]))
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g)
#' extract_dn(caim, sky_points)
#'
#' # ImageJ can be used to digitize points.
#' # See calc_zenith_raster_coord() for details.
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' plot(caim)
#' path <- system.file("external/points_over_perimeter.csv",
#'                      package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[,c(ncol(sky_points)-1, ncol(sky_points))]
#' colnames(sky_points) <- c("col", "row")
#' head(sky_points)
#' v <- cellFromRowCol(caim, sky_points$row, sky_points$col) %>%
#'   xyFromCell(caim, .) %>% vect()
#' plot(v, add = TRUE, col = 2)
#' extract_dn(caim, sky_points)
extract_dn <- function(r, sky_points, use_window = TRUE) {
  stopifnot(is.data.frame(sky_points))
  stopifnot(is(r, "SpatRaster"))

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  if (use_window) {
    xy <-  terra::xyFromCell(r, cells)
    r_smooth <- terra::focal(r, 3, "mean")
    dn <- terra::extract(r_smooth, xy)[,]
    if(any(is.na(sky_points))) {
      warning("The kernel created NA values near the horizon.")
    }
  } else {
    dn <- terra::extract(r, xy)[,]
  }
  .names <- c(colnames(sky_points), names(r))
  sky_points <- cbind(sky_points, dn)
  colnames(sky_points) <- .names
  sky_points
}
