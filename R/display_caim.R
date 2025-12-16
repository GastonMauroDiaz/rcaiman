#' Display a canopy image
#'
#' Wrapper for [EBImage::display()] that streamlines the visualization of
#' canopy images, optionally overlaying binary masks and segmentation borders.
#' It is intended for quick inspection of processed or intermediate results in
#' a graphical viewer.
#'
#' @param caim [terra::SpatRaster-class]. Typically the output of [read_caim()].
#'   Can be multi- or single-layer.
#' @param seg Segmentation map typically created with functions such as
#'   [equalarea_segmentation()], [skygrid_segmentation()], [ring_segmentation()]
#'   or [sector_segmentation()], but any raster with integer segment labels is
#'   accepted.
#' @param sun_row_col numeric `data.frame` with the estimated sun‑disk
#'   position in image coordinates. See [row_col_from_zenith_azimuth()].
#' @param sun_disk_size numeric vector of length one. Sun disk size in pixels.
#'
#' @inheritParams compute_canopy_openness
#' @inheritParams extract_dn
#'
#' @return Invisible `NULL`. Called for side effects (image viewer popup).
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
#' x11()
#' # plot(caim$Blue)
#' # sun_angles <- click(c(z, a), 1) %>% as.numeric()
#' sun_angles <- c(z = 49.5, a = 27.4) #taken with above lines then hardcoded
#'
#' sun_row_col <- row_col_from_zenith_azimuth(z, a,
#'                                            sun_angles["z"],
#'                                            sun_angles["a"])
#'
#' bin <- binarize_with_thr(r$Blue, thr_isodata(r$Blue[]))
#' seg <- equalarea_segmentation(z, a, 200)
#'
#' # color image
#' display_caim(caim)
#' # greyscale
#' display_caim(caim$Blue)
#' # to check binarization quality (press `h` with the pointer on the image to
#' # learn useful hotkeys)
#' display_caim(caim$Blue, bin)
#' # to see the segments
#' display_caim(seg = seg)
#' # to see the marks
#' display_caim(caim, sampling_points = sky_points, sun_row_col = sun_row_col)
#' }
display_caim <- function(caim = NULL,
                         bin = NULL,
                         seg = NULL,
                         sampling_points = NULL,
                         sun_row_col = NULL,
                         sun_disk_size = 9) {

  .check_vector(sun_disk_size, "numeric", 1, sign = "positive")

  .this_requires_EBImage()

  if (!is.null(caim)) {
    .assert_spatraster(caim)
  }
  if (!is.null(bin)) {
    .assert_logical_mask(bin)
  }
  if (!is.null(seg)) {
    .assert_single_layer(seg)
    # Highlight region borders
    if (!is.null(seg)) {
      laplacian <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), nrow = 3)
      seg <- terra::focal(seg, laplacian)
      seg <- seg != 0
    }
  }
  if (!is.null(sampling_points)) {
    .check_sky_points(sampling_points)
    # create raster
    sampling_points <- cbind(sampling_points, dn = 1)
    sampling_points <- interpolate_planar(sampling_points, caim[[1]], k = 1, p = 1, rmax = 1.5, col_id = 3)
    sampling_points <- is.na(sampling_points)
  }
  if (!is.null(sun_row_col)) {
    .check_sky_points(sun_row_col)
    # create raster
    if (any(sun_row_col < 1)) {
      warning(paste0("The ", names(sun_row_col)[which(sun_row_col < 1)]," number of `sun_row_col` was forced to 1."))
      sun_row_col[sun_row_col < 1] <- 1
    }
    sun_row_col <- cbind(sun_row_col, dn = 1)
    sun_row_col <- interpolate_planar(sun_row_col, caim[[1]], k = 1, p = 1, rmax = sun_disk_size, col_id = 3)
    sun_row_col <- is.na(sun_row_col)
  }

  if (!is.null(caim) && terra::nlyr(caim) <= 3 && is.null(bin) &&
      is.null(seg) && (!is.null(sampling_points) | !is.null(sun_row_col))) {
    caim <- normalize_minmax(caim)
    if (!is.null(sun_row_col)) {
      caim <- paint_with_mask(caim, sun_row_col, color = "yellow")
    }
    if (!is.null(sampling_points)) {
      caim <- paint_with_mask(caim, sampling_points, color = "red")
    }
    x <- EBImage::Image(caim, dim(caim), colormode = "color")
  } else {
    layers <- list()

    if (!is.null(caim)) {
      layers <- c(layers, normalize_minmax(caim))
    }
    if (!is.null(bin)) {
      layers <- c(layers, bin)
    }
    if (!is.null(seg)) {
      layers <- c(layers, seg)
    }

    if (length(layers) > 0) {
      if (is.null(bin) && is.null(seg)) {
        if (length(layers) == 1)
          if (terra::nlyr(caim) == 3) {
            x <- EBImage::Image(rast(layers), dim(caim), colormode = "color")
          } else {
            x <- EBImage::Image(rast(layers), dim(caim))
          }
      } else {
        x <- do.call(c, layers)
        x <- EBImage::Image(x, dim(x))
      }
    } else {
      warning("Nothing to display")
      x <- NULL
    }
  }
  EBImage::display(x)
}

