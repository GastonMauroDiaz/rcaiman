#' Read a canopy image from a file
#'
#' Wrapper function for \code{\link[raster]{raster}}.
#'
#'
#' Run \code{read_caim()} to obtain an example of a hemispherical photo taken in
#' non-diffuse light conditions in a \emph{Nothofagus pumilio} forest from
#' Argentina with a FC-E9 auxiliary lens attached to a Nikon Coolpix 5700.
#'
#' Since this function aims to read born-digital color photographs, RGB-JPEG and
#' RGB-TIFF are expected as input. To read a region of the file use
#' \code{upper_left}, \code{width}, and \code{height}. The \code{upper_left}
#' parameter indicates the pixels coordinates of the upper left corner of the
#' region of interest (ROI). Those coordinates should be in the raster
#' coordinates system, which works like a spreadsheet, i.e, when you go down
#' through the vertical axis, the "row" number increase (*IMPORTANT: column and
#' row must be provided instead of row and column*). The \code{width}, and
#' \code{height} parameters indicate the size of the boxy ROI. I recommend using
#' \href{https://imagej.nih.gov/ij/}{ImageJ} to obtain this parameters, but any
#' image editor can be used, such as GIMP and Adobe Photoshop.
#'
#'
#' @param path_to_file Character vector of length one. Path to a JPEG or TIFF
#'   file. The function will return a data example (see details) if no arguments
#'   are provided.
#' @param upper_left An integer vector of length two (see details).
#' @param width,height An integer vector of length one (see details).
#'
#' @return \code{\linkS4class{RasterBrick}}.
#' @export
#'
#' @seealso write_caim
#'
#' @examples
#' # This is the example image
#' r <- read_caim()
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
setGeneric("read_caim", function(path_to_file, upper_left = NULL, width = NULL,
                            height = NULL) {
  standardGeneric("read_caim")
})

#' @describeIn read_caim Provide the path to a file. If The file is stored in
#'   the working directory, just provide the file name. File extension should be
#'   included in the file name.
setMethod(
  "read_caim",
  signature(path_to_file = "character"),
  function(path_to_file, upper_left, width, height) {
    file_ext <- toupper(extension(path_to_file))

    if (file_ext %in% c(".JPG", ".JPEG", ".TIF", ".TIFF")) {
      r <- brick(path_to_file)
      if (nlayers(r) != 3) stop("The photograph should have three layers.")
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

      xmn <- xFromCol(r, upper_left[1])
      xmx <- xFromCol(r, upper_left[1] + width)
      ymx <- yFromRow(r, upper_left[2])
      ymn <- yFromRow(r, upper_left[2] + height)

      if (any(is.na(xmn), is.na(xmx), is.na(ymn), is.na(ymx))) {
        stop(
          paste(
            "The selection is outside the picture border, please",
            "review upper_left, height, and width."
          )
        )
      }

      e <- extent(xmn, xmx, ymn, ymx)
      r <- crop(r, e)
      extent(r) <- extent(0, ncol(r), 0, nrow(r))
      projection(r) <- NA
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


#' Write canopy image
#'
#' Wrapper function for \code{\link[raster]{writeRaster}}.
#'
#' @param caim \linkS4class{Raster}.
#' @param path Character vector of length one. Path for writing the image.
#' @param bit_depth Numeric vector of length one.
#'
#' @export
#'
#' @seealso write_bin
#'
#' @examples
#' \dontrun{
#' require(magrittr)
#' caim <- read_caim() %>% normalize(., 0, 255)
#' write_caim(caim * 2^8, "test_8bit", 8)
#' write_caim(caim * 2^16, "test_16bit", 16)
#' }
write_caim <- function(caim, path, bit_depth) {

  if (!any(bit_depth == 16, bit_depth == 8)) {
    stop("bit_depth should be 8 or 16.")
  }

  projection(caim) <- NA
  extent(caim) <- extent(0, ncol(caim), 0, nrow(caim))

  file_name <- basename(path)
  extension(file_name) <- "tif"

  if (bit_depth == 8) {
    suppressWarnings(
      writeRaster(caim, file.path(dirname(path), file_name),
                format = "GTiff", datatype = "INT1U", overwrite = TRUE)
    )
  } else {
    suppressWarnings(
      writeRaster(caim, file.path(dirname(path), file_name),
                  format = "GTiff", datatype = "INT2U", overwrite = TRUE)
    )
  }
}





