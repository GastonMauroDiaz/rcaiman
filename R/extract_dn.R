#' Extract digital numbers
#'
#' It is a wrapper function around \code{\link[terra]{extract}}.
#'
#'
#' @inheritParams fisheye_to_equidistant
#' @param img_points The result of a call to \code{\link{extract_sky_points}},
#'   or an object of the same class and structure.
#' @inheritParams extract_rl
#' @inheritParams extract_feature
#'
#' @return An object of the class \emph{data.frame}, which is the argument
#'   \code{img_points} with an added column per each layer from \code{r}. The
#'   layer names are used to name the new columns. If a function is provided as
#'   the \code{fun} argument, the result will be summarized per column using the
#'   provided function, and the \emph{row} and \emph{col} information will be
#'   omitted. Moreover, if \code{r} is an RGB image, a \linkS4class{color} will
#'   be returned instead of a \emph{data.frame}. The latter feature is useful
#'   for obtaining  the \code{sky_blue} argument for \code{\link{enhance_caim}}.
#' @export
#'
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- gbc(caim$Blue)
#' bin <- apply_thr(r, thr_isodata(r[]))
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_dn(caim, sky_points)
#' head(sky_points)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points
#' }
#'
#' # ImageJ can be used to digitize points.
#' # See calc_zenith_raster_coord() for details.
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' plot(caim)
#' path <- system.file("external/points_over_perimeter.csv",
#'                      package = "rcaiman")
#' img_points <- read.csv(path)
#' img_points <- img_points[,c(ncol(img_points), ncol(img_points)-1)]
#' colnames(img_points) <- c("row", "col")
#' head(img_points)
#' v <- cellFromRowCol(caim, img_points$row, img_points$col) %>%
#'   xyFromCell(caim, .) %>% vect()
#' plot(v, add = TRUE, col = 2)
#' extract_dn(caim, img_points, fun = median)
extract_dn <- function(r, img_points, use_window = TRUE, fun = NULL) {
  stopifnot(is.data.frame(img_points))
  stopifnot(is(r, "SpatRaster"))

  cells <- terra::cellFromRowCol(r, img_points$row, img_points$col)
  if (use_window) {
    xy <-  terra::xyFromCell(r, cells)
    r_smooth <- terra::focal(r, 3, "mean")
    dn <- terra::extract(r_smooth, xy)[,]
    if(any(is.na(img_points))) {
      warning("The kernel produced NA values near the horizon.")
    }
  } else {
    dn <- terra::extract(r, xy)[,]
  }
  .names <- c(colnames(img_points), names(r))
  img_points <- cbind(img_points, dn)
  colnames(img_points) <- .names
  if (!is.null(fun)) {
      if (ncol(img_points) > 3) {
        img_points <- apply(img_points[,-(1:2)], 2, fun)
      } else {
        img_points <- fun(img_points[,3])
      }

      if (all(names(r) == names(read_caim()))) {
        img_points <- colorspace::sRGB(img_points[1],
                                       img_points[2],
                                       img_points[3])
      }
    }
  img_points
}
