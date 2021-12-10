#' Expand non-circular
#'
#' Expand a non-circular hemispherical photograph.
#'
#' @param caim \linkS4class{RasterBrick}. The return of a call to
#'   \code{\link{read_caim}}.
#' @inheritParams azimuth_image
#' @param zenith_colrow Numeric vector of length two. Raster coordinates of the
#'   zenith. See \code{\link{calc_zenith_raster_coordinates}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    my_file <- path.expand("~/DSC_2881.JPG")
#'    download.file("https://osf.io/x8urg/download", my_file,
#'                method = "auto", mode = "wb"
#'    )
#'
#'    r <- read_caim(file.path(path, "DSC_2881.JPG"))
#'    diameter <- calc_diameter(lens("Nikkor_10.5_mm"), 1202, 53)
#'    zenith_colrow <- c(1503, 998)
#'    z <- zenith_image(diameter, lens("Nikkor_10.5_mm"))
#'    r <- expand_noncircular(r, z, zenith_colrow)
#' }
expand_noncircular <-  function (caim, z, zenith_colrow) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(zenith_colrow) == "numeric")
  stopifnot(length(zenith_colrow) == 2)

  zenith_xy <- zenith_colrow
  zenith_xy[2] <- nrow(caim) - zenith_colrow[2]

  # locate the center of caim over the center of z
  center <- ncol(z) / 2
  xmn <- center - (ncol(caim) / 2)
  xmx <- center + (ncol(caim) / 2)
  ymn <- center - (nrow(caim) / 2)
  ymx <- center + (nrow(caim) / 2)
  e <- extent(xmn, xmx, ymn, ymx)
  extent(caim) <- e

  # shift the center of caim according with zenith_xy
  delta_x <- zenith_xy[1] - ncol(caim) / 2
  delta_y <- zenith_xy[2] - nrow(caim) / 2
  xmn <- xmin(caim) + delta_x
  xmx <- xmax(caim) + delta_x
  ymn <- ymin(caim) + delta_y
  ymx <- ymax(caim) + delta_y
  extent(caim) <- extent(xmn, xmx, ymn, ymx)

  r <- extend(caim, z, value = NA)
  extent(r) <- extent(0, ncol(r), 0, nrow(r))
  r

}
