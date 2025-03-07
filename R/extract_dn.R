#' Extract digital numbers
#'
#' Wrapper function around [terra::extract()].
#'
#'
#' @inheritParams fisheye_to_equidistant
#' @param img_points The output of [extract_sky_points()],
#'   or an object of the same class and structure.
#' @inheritParams extract_rl
#' @inheritParams extract_feature
#'
#' @return An object of the class *data.frame*. It is the argument
#'   `img_points` with an added column per each layer from `r`. The
#'   layer names are used to name the new columns. If a function is provided as
#'   the `fun` argument, the result will be summarized per column using the
#'   provided function, and the *row* and *col* information will be
#'   omitted. Moreover, if `r` is an RGB image, a [color-class] will
#'   be returned instead of a *data.frame*. The latter feature is useful
#'   for obtaining  the `sky_blue` argument for [enhance_caim()].
#' @export
#'
#' @note
#' The [point selection tool of ‘ImageJ’
#' software](https://imagej.net/ij/docs/guide/146-19.html#sec:Multi-point-Tool)
#' can be used to manually digitize points and create a CSV file from which to read
#' coordinates (see Examples). After digitizing the points on the image, use the
#' dropdown menu Analyze>Measure to open the Results window. To obtain the CSV
#' file, use File>Save As...
#'
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' bin <- apply_thr(r, thr_isodata(r[]))
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_dn(caim, sky_points)
#' head(sky_points)
#'
#' # See fit_cie_sky_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' extract_dn(caim, sky_points, fun = median)
#' }
extract_dn <- function(r, img_points, use_window = TRUE, fun = NULL) {
  stopifnot(is.data.frame(img_points))
  stopifnot(ncol(img_points) == 2)
  stopifnot(is(r, "SpatRaster"))

  cells <- terra::cellFromRowCol(r, img_points$row, img_points$col)
  xy <-  terra::xyFromCell(r, cells)
  if (use_window) {
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
