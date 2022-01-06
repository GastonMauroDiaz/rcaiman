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
#' @family Lens functions
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
#'    r <- read_caim(my_file)
#'    diameter <- calc_diameter(lens("Nikkor_10.5_mm"), 1202, 53)
#'    zenith_colrow <- c(1503, 998)
#'    z <- zenith_image(diameter, lens("Nikkor_10.5_mm"))
#'    r <- expand_noncircular(r, z, zenith_colrow)
#'    plot(r)
#' }
expand_noncircular <-  function (caim, z, zenith_colrow) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(zenith_colrow) == "numeric")
  stopifnot(length(zenith_colrow) == 2)

  zenith_xy <- c(zenith_colrow[1], nrow(caim) - zenith_colrow[2])
  delta_x <-  (ncol(caim) / 2) - zenith_xy[1]
  delta_y <-  (nrow(caim) / 2) - zenith_xy[2]
  center <- ncol(z) / 2
  xmn <- center - (ncol(caim) / 2) - delta_x
  xmx <- center + (ncol(caim) / 2) - delta_x
  ymn <- center - (nrow(caim) / 2) + delta_y
  ymx <- center + (nrow(caim) / 2) + delta_y
  e <- extent(xmn, xmx, ymn, ymx)
  extent(caim) <- e

  ze <- extent(z) * 1.5
  r <- extend(caim, z, value = NA)
  extent(r) <- alignExtent(extent(r), z)
  r <- crop(r, z)
  extent(r) <- extent(0, ncol(r), 0, nrow(r))
  r
}
