#' Read a canopy image from a file
#'
#' Wrapper function for [terra::rast()].
#'
#' Run `read_caim()` to obtain an example of a hemispherical photo taken in
#' non-diffuse light conditions in a *Nothofagus pumilio* forest with a FC-E9
#' auxiliary lens attached to a Nikon Coolpix 5700.
#'
#' Since this function aims to read born-digital color photographs, RGB-JPEG and
#' RGB-TIFF are the expected input. However, since this function is a wrapper
#' for [terra::rast()], format compatibility is heritages from it.
#'
#' Use `upper_left`, `width`, and `height` to read a particular region from the
#' file. Although any image editor can be used to obtain those parameters, this
#' function was tested with [‘ImageJ’](https://imagej.nih.gov/ij/), so it might
#' be wise to use it.
#'
#' **TIP**: For obtaining `upper_left`, `width`, and
#' `height`, open the image on the Fiji distro of ImageJ, draw a rectangular
#' selection, and go to Edit>Selection>Specify. The same workflow may work with
#' other distros.
#'
#' @param path Character vector of length one. Path to an image, including file
#'   extension. The function will return a data example if no arguments are
#'   provided.
#' @param upper_left An integer vector of length two. The pixels coordinates of
#'   the upper left corner of a region of interest (ROI). These coordinates
#'   should be in the raster coordinates system. This system works like a
#'   spreadsheet, i.e, when going down through the vertical axis, the *row*
#'   number increases (**IMPORTANT**: column and row must be provided instead of
#'   row and column, as is the norm for objects from the class *data.frame* and
#'   others alike)
#' @param width,height An integer vector of length one. The size of the boxy ROI
#'   whose upper left corner is the `upper_left` argument.
#'
#' @return An object from class [SpatRaster-class] with its layers named
#'   *Red*, *Green*, and *Blue* when a born-digital color
#'   photographs is provided as input.
#'
#' @note
#'
#' The example image was obtained with this code:
#'
#' ````
#' zenith_colrow <- c(1290, 988)
#' z <- zenith_image(745*2, lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- read_caim_raw("DSCN4606.NEF", z, a, zenith_colrow, radius = 300)
#' z <- zenith_image(ncol(r), lens())
#' r <- correct_vignetting(r, z, c(0.0638, -0.101))
#' r <- c(mean(r$Y, r$M), r$G, r$C)
#' r <- normalize(r, -1)
#' write_caim(r*2^16-2, "example.tif", 16)
#' ````
#'
#' @export
#' @family Tool Functions
#'
#'
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' #plot(caim$Blue)
read_caim <- function(path = NULL, upper_left = NULL, width = NULL,
                      height = NULL) {

  if (is.null(path)) {
    path <- system.file("external/example.tif", package = "rcaiman")
  }

  suppressWarnings(r <- terra::rast(path))
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  # https://spatialreference.org/ref/sr-org/7589/
  terra::crs(r) <- "epsg:7589"


  if (all(!is.null(upper_left), !is.null(height), !is.null(width))) {
    if (length(upper_left) != 2) {
      stop("upper_left should be a numeric vector of length two")
    }
    if (any(upper_left == c(0, 0))) {
      stop("upper_left should be c(1, 1) instead of c(0,0).")
    }
    if (all(length(height) != 1, as.integer(height) == height)) {
      stop("height should be a one-lenght integer")
    }
    if (all(length(width) != 1, as.integer(width) == width)) {
      stop("width should be a one-lenght integer")
    }

    xmn <- terra::xFromCol(r, upper_left[1])
    xmx <- terra::xFromCol(r, upper_left[1] + width)
    ymx <- terra::yFromRow(r, upper_left[2])
    ymn <- terra::yFromRow(r, upper_left[2] + height)

    if (any(is.na(xmn), is.na(xmx), is.na(ymn), is.na(ymx))) {
      stop(
        paste(
          "The selection is outside the picture border, please",
          "review upper_left, height, and width."
        )
      )
    }
    e <- terra::ext(xmn, xmx, ymn, ymx)
    r <- terra::crop(r, e)
    terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  }

  file_ext <- filenamer::as.filename(path)$ext %>% toupper()

  if (file_ext %in% c("JPG", "JPEG", "TIF", "TIFF")) {
    if (terra::nlyr(r) == 3) names(r) <- c("Red", "Green", "Blue")
  } else {
    if (terra::nlyr(r) == 3) {
      names(r) <- c("Red", "Green", "Blue")
      warning("Layers were named presuming an RGB file")
    }
  }
  r
}
