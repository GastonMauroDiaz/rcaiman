#' Obtain zenith and azimuth angles from row and col numbers
#'
#' @note
#' Use the `lens_coef` argument to calculate coordinates below the horizon.
#'
#' @inheritParams ootb_mblt
#' @inheritParams zenith_image
#' @param row_col Numeric vector of length two. Row and col numbers.
#'
#' @export
#'
#' @family HSP Functions
#' @return Numeric vector of length two.
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' zenith_azimuth_from_row_col(z, a, c(501, 750))
zenith_azimuth_from_row_col <- function(z, a, row_col, lens_coef = NULL) {

  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(z, a)
  stopifnot(is.numeric(row_col))
  stopifnot(length(row_col) == 2)

  v <- cellFromRowCol(z, row_col[1], row_col[2]) %>%
    xyFromCell(z, .) %>% vect()

  if (is.null(lens_coef)) {
    l <- list(c(terra::extract(z, v)[1,2],
                terra::extract(a, v)[1,2]), row_col)
  } else {
    z <- terra::deepcopy(z)

    #get azimuth
    e <- terra::ext(z)
    terra::ext(z) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
    xy <- terra::cellFromRowCol(z, row_col[1], row_col[2]) %>%
      terra::xyFromCell(z, .)
    tr <- pracma::cart2pol(as.numeric(xy))
    #get relative radius
    rr <- tr[2] * 2/pi
    #invert
    zs <- seq(0,150, 0.1)
    rrs <- calc_relative_radius(zs, lens_coef)
    z_from_rr <- suppressWarnings(splinefun(rrs, zs))
    zenith <- z_from_rr(rr)
    l <- list(c(zenith, terra::extract(a, v)[1,2]), row_col)
  }
  names(l) <- c("zenith_azimuth", "row_col")
  l
}
