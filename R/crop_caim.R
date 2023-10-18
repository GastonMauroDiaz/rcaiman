#' Crop a canopy image from a file
#'
#' Function that complements [read_caim()] and [read_caim_raw()]
#'
#' @param r [SpatRaster-class]
#' @inheritParams read_caim
#'
#' @return [SpatRaster-class]
#' @export
#' @examples
#' caim <- read_caim()
#' ncell(caim)
#' caim <- crop_caim(caim, c(231,334), 15, 10)
#' ncell(caim)
crop_caim <- function(r, upper_left = NULL, width = NULL,
                      height = NULL) {

  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  # https://spatialreference.org/ref/sr-org/7589/
  terra::crs(r) <- "epsg:7589"

  # START code from read_caim
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
  # END code from read_caim()
  r
}
