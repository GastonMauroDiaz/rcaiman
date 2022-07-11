#' Expand non-circular
#'
#' Expand a non-circular hemispherical photograph.
#'
#' @param caim \linkS4class{SpatRaster}. The return of a call to
#'   \code{\link{read_caim}}.
#' @inheritParams azimuth_image
#' @param zenith_colrow Numeric vector of length two. Raster coordinates of the
#'   zenith. See \code{\link{calc_zenith_raster_coordinates}}.
#'
#' @family Lens functions
#'
#' @return An object of class \linkS4class{SpatRaster} that is the result of
#'   adding margins (\code{NA} pixels) to \code{caim}. The zenith point depicted
#'   in the picture should be in the center of the image or very close to it.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    my_file <- file.path(tempdir(), "DSC_2881.JPG")
#'    download.file("https://osf.io/x8urg/download", my_file,
#'                method = "auto", mode = "wb"
#'    )
#'
#'    r <- read_caim(my_file)
#'    diameter <- calc_diameter(lens("Nikkor_10.5_mm"), 1202, 53)
#'    zenith_colrow <- c(1503, 998)
#'    z <- zenith_image(diameter, lens("Nikkor_10.5_mm"))
#'    r <- expand_noncircular(r, z, zenith_colrow)
#'    plot(r)
#'    plot(is.na(r$Red), add = TRUE, alpha = 0.5)
#' }
expand_noncircular <-  function (caim, z, zenith_colrow) {
  .is_single_layer_raster(z, "z")
  stopifnot(class(zenith_colrow) == "numeric")
  stopifnot(length(zenith_colrow) == 2)

  zenith_xy <- c(zenith_colrow[1], nrow(caim) - zenith_colrow[2])
  delta_x <-  (ncol(caim) / 2) - zenith_xy[1]
  delta_y <-  (nrow(caim) / 2) - zenith_xy[2]
  center <- ncol(z) / 2
  if (ncol(caim) > nrow(caim)) {
    xmn <- center - (ncol(caim)/2) - delta_x
    xmx <- center + (ncol(caim)/2) - delta_x
  } else {
    xmn <- center - (ncol(caim)/2) + delta_x
    xmx <- center + (ncol(caim)/2) + delta_x
  }
  ymn <- center - (nrow(caim) / 2) + delta_y
  ymx <- center + (nrow(caim) / 2) + delta_y
  e <- terra::ext(xmn, xmx, ymn, ymx)
  terra::ext(caim) <- e

  ze <- terra::ext(z) * 1.5
  r <- terra::extend(caim, z)
  terra::ext(r) <- terra::align(terra::ext(r), z)
  r <- terra::crop(r, z)
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  # names(r) <- "Expanded non-circular image"
  r
}
