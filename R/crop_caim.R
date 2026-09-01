#' Crop a canopy image
#'
#' Extracts a rectangular region of interest (ROI) from a canopy image. This
#' function complements [read_caim()].
#'
#' @param r [terra::SpatRaster-class].
#' @param upper_left numeric vector of length two. Pixel coordinates of the
#'   upper-left corner of the ROI, in the format `c(column, row)`.
#' @param width,height numeric vector of length one. Size (in pixels) of the
#'   rectangular ROI to read.
#'
#'
#' @section Selecting a Region of Interest:
#' Arguments `upper_left`, `width`, and `height` are expressed in raster
#' coordinates, similar to a spreadsheet layout: **columns first, then rows**.
#' In other words, specify coordinates as `c(column, row)`, **not** `c(row,
#' column)`, which is typical in `data.frame` objects.
#'
#' While any image editor can be used to obtain these values, this function was
#' tested with [ImageJ](https://imagej.net/ij/), particularly the Fiji
#' distribution. A recommended workflow:
#' 1. Open the image in Fiji.
#' 2. Draw a rectangular selection.
#' 3. Go to *Edit > Selection > Specify...* to read `upper_left`, `width`, and `height`.
#'
#' @return [terra::SpatRaster-class] object containing the same layers and values
#'   as `r` but restricted to the selected ROI, preserving all other properties.
#'
#' @note `rcaiman` uses terra without geographic semantics: rasters are kept with
#'   unit resolution (cell size = 1) and a standardized extent
#'   `ext(0, ncol, 0, nrow)` with CRS EPSG:7589 (see [set_rcaiman_geometry()]).
#'
#' @export
#'
#' @examples
#' caim <- read_caim()
#' ncell(caim)
#' caim <- crop_caim(caim, c(231,334), 15, 10)
#' ncell(caim)
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_jpeg.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry, "exif_jpeg", "delta_t")
#'
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path)
#' zenith_col_row <- spec$zenith_col_row
#' horizon_radius <- spec$horizon_radius
#' upper_left <- zenith_col_row - horizon_radius
#' side <- horizon_radius * 2
#' plot(crop_caim(caim, upper_left, side, side))
#' }
crop_caim <- function(r, upper_left = NULL, width = NULL,
                      height = NULL) {

  .assert_spatraster(r)
  .check_vector(upper_left, "numeric", 2, allow_null = TRUE, sign = "positive")
  .check_vector(width, "numeric", 1, allow_null = TRUE, sign = "positive")
  .check_vector(height, "numeric", 1, allow_null = TRUE, sign = "positive")

  r <- set_rcaiman_geometry(r)

  if (all(!is.null(upper_left), !is.null(height), !is.null(width))) {
    xmn <- terra::xFromCol(r, upper_left[1])
    xmx <- terra::xFromCol(r, upper_left[1] + width)
    ymx <- terra::yFromRow(r, upper_left[2])
    ymn <- terra::yFromRow(r, upper_left[2] + height)

    if (any(is.na(xmn), is.na(xmx), is.na(ymn), is.na(ymx))) {
      stop(
        "The selection is outside the picture border, review `upper_left`, `height`, and `width`."
      )
    }
    e <- terra::ext(xmn, xmx, ymn, ymx)
    r <- terra::crop(r, e)
    r <- set_rcaiman_geometry(r)
  }

  r
}
