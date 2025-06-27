#' Extract digital numbers
#'
#' Wrapper function around [terra::extract()].
#'
#'
#' @param r [SpatRaster-class]. The image from which the argument `sky_points`
#'   was obtained.
#' @param sky_points The output of [extract_sky_points()]
#'   or an object of the same class and structure.
#' @param use_window Logical vector of length one. If `TRUE`, a window of \eqn{3
#'   \times 3} pixels will be used to extract the digital number from `r`.
#' @inheritParams extract_feature
#'
#' @return An object of the class *data.frame*. It is the argument
#'   `sky_points` with an added column per each layer from `r`. The
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
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # See fit_cie_sky_model() for details on below file
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
#' }
extract_dn <- function(r, sky_points, use_window = TRUE, fun = NULL) {
  stopifnot(is.data.frame(sky_points))
  stopifnot(ncol(sky_points) == 2)
  stopifnot(is(r, "SpatRaster"))

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  xy <-  terra::xyFromCell(r, cells)
  if (use_window) {
    r_smooth <- terra::focal(r, 3, "mean")
    dn <- terra::extract(r_smooth, xy)[,]
    if(any(is.na(sky_points))) {
      warning("The kernel produced NA values near the horizon.")
    }
  } else {
    dn <- terra::extract(r, xy)[,]
  }
  sp_names <- c(colnames(sky_points), names(r))
  sky_points <- cbind(sky_points, dn)
  colnames(sky_points) <- sp_names
  if (!is.null(fun)) {
      if (ncol(sky_points) > 3) {
        sky_points <- apply(sky_points[,-(1:2)], 2, fun)
      } else {
        sky_points <- fun(sky_points[,3])
      }

      if (terra::nlyr(r) == 3) {
        if (all(names(r) == names(read_caim()))) {
          sky_points <- colorspace::sRGB(sky_points[1],
                                         sky_points[2],
                                         sky_points[3])
        }
      }
    }
  sky_points
}
