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
#' # See fit_cie_sky_model() for details on below file
#' path <- system.file("external/sky_points.gpkg",
#'                     package = "rcaiman")
#' sky_points <- terra::vect(path)
#' sky_points <- terra::extract(caim, sky_points, cells = TRUE)
#' sky_points <- terra::rowColFromCell(caim, sky_points$cell) %>%
#'   as.data.frame()
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 15), "thr_isodata")
#' bin <- bin & mask_hs(z, 0, 85)
#' mx <- optim_normalize(caim, bin)
#' caim <- normalize(caim, mx = mx, force_range = TRUE)
#' m <- !is.na(z)
#'
#' sky_blue <- extract_dn(caim, sky_points, fun = median)
#' as(sky_blue, "polarLAB")
#' ecaim <- enhance_caim(caim, m, sky_blue = sky_blue)
#' plot(ecaim)
#' .refine_sky_blue <- function(chroma) {
#'   ecaim <- enhance_caim(caim, m, polarLAB(50, chroma, 293))
#'   z <- -extract_dn(ecaim, sky_points, fun = sum)
#'   if (is.na(z)) z <- 0
#'   z
#' }
#' chroma <- seq(0, 1, 0.05) * 100
#' z <- Map(.refine_sky_blue, chroma) %>% unlist()
#' plot(chroma, z, type = "l")
#'
#' fit <- bbmle::mle2(.refine_sky_blue, list(chroma = 17),
#'                    method = "L-BFGS-B", lower = 0, upper = 100)
#' fit@details$convergence
#' summary(fit)
#' fit@coef
#' ecaim <- enhance_caim(caim, m, polarLAB(50, 17, 293))
#' plot(ecaim)
#' apply_thr(ecaim, thr_isodata(ecaim[m])) %>% plot()
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
