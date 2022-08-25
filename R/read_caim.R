#' Read a canopy image from a file
#'
#' Wrapper function for \code{\link[terra]{rast}}.
#'
#'
#' Run \code{read_caim()} to obtain an example of a hemispherical photo taken in
#' non-diffuse light conditions in a \emph{Nothofagus pumilio} forest with a
#' FC-E9 auxiliary lens attached to a Nikon Coolpix 5700.
#'
#' Since this function aims to read born-digital color photographs, RGB-JPEG and
#' RGB-TIFF are expected as input. Use \code{upper_left}, \code{width}, and
#' \code{height} to read a region of the file. The \code{upper_left} parameter
#' indicates the pixels coordinates of the upper left corner of the region of
#' interest (ROI). These coordinates should be in the raster coordinates system,
#' which works like a spreadsheet, i.e, when you go down through the vertical
#' axis, the \emph{row} number increases (\strong{IMPORTANT: column and row must
#' be provided instead of row and column as in objects from the class data.frame
#' and others alike}). The \code{width} and \code{height} parameters indicate
#' the size of the boxy ROI. I recommend using
#' \href{https://imagej.nih.gov/ij/}{‘ImageJ’} to obtain these parameters, but
#' any image editor can be used, such as ‘GIMP’ or ‘Adobe Photoshop’.
#'
#'
#' @param path_to_file Character vector of length one. Path to a JPEG or TIFF
#'   file. The function will return a data example if no arguments
#'   are provided.
#' @param upper_left An integer vector of length two.
#' @param width,height An integer vector of length one.
#'
#' @return An object from class \linkS4class{SpatRaster} with its layers named
#'   \emph{Red}, \emph{Green}, and \emph{Blue}.
#' @export
#' @family Tools Functions
#'
#'
#' @examples
#' # This is the example image
#' r <- read_caim()
#' plotRGB(r)
#'
#' # This is also the example
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' # the zenith raster coordinates can be easily transformed to the "upper_left"
#' # argument by subtracting from it the radius expressed in pixels.
#' zenith_colrow <- c(1280, 960)
#' diameter_px <- 1490
#' r <- read_caim(path,
#'                upper_left = zenith_colrow - diameter_px/2,
#'                width = diameter_px,
#'                height = diameter_px)
#' plotRGB(r)
setGeneric("read_caim", function(path_to_file, upper_left = NULL, width = NULL,
                            height = NULL) {
  standardGeneric("read_caim")
})

#' @describeIn read_caim Provide the path to a file. If the file is stored in
#'   the working directory, providing the file name is enough. File extension
#'   should be included in the file name.
setMethod(
  "read_caim",
  signature(path_to_file = "character"),
  function(path_to_file, upper_left, width, height) {
    file_ext <- filenamer::as.filename(path_to_file)$ext %>% toupper()

    if (file_ext %in% c("JPG", "JPEG", "TIF", "TIFF")) {
      suppressWarnings(r <- terra::rast(path_to_file))
      if (terra::nlyr(r) != 3) stop("The photograph should have three layers.")
      terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
      # https://spatialreference.org/ref/sr-org/7589/
      terra::crs(r) <- "epsg:7589"
    } else {
      stop("The file extension does not correspond to a JPEG or TIFF file.")
    }

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

    names(r) <- c("Red", "Green", "Blue")

    r
  }
)

#' @describeIn read_caim It returns an example (see details).
setMethod(
  "read_caim",
  signature(path_to_file = "missing"),
  function(path_to_file) {
    path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
    read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
  }
)

