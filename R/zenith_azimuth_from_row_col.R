#' Obtain zenith and azimuth angles from row and col numbers
#'
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams zenith_image
#' @param row,col Numeric vector. Row or column numbers.
#'
#' @export
#'
#' @return An object of the class _data.frame_.
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' zenith_azimuth_from_row_col(z, a, 499, 751)
zenith_azimuth_from_row_col <- function(z, a, row, col) {

  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(z, a)
  stopifnot(length(row) == length(col))

  i <- terra::cellFromRowCol(z, row, col)
  u <- is.na(z[i])[,]

  if (any(u)) {
    size <- 100
    sky_points <- expand.grid(row = seq(1, nrow(z), length.out = size) %>%
                                round(),
                              col = seq(1, nrow(z), length.out = size) %>%
                                round())
    sky_points <- extract_dn(z, sky_points, use_window = FALSE)
    sky_points <- sky_points[!is.na(sky_points[,3]), c("row", "col")]
    sky_points <- extract_rel_radiance(z, z, a, sky_points,
                                       use_window = FALSE)$sky_points

    fit <- spatial::surf.ls(x = sky_points[, "row"],
                            y = sky_points[, "col"],
                            z = sky_points[, "z"],
                            np = 3)
    zenith <- predict(fit, row, col) %>% round()

    fit <- spatial::surf.ls(x = sky_points[, "row"],
                            y = sky_points[, "col"],
                            z = sky_points[, "a"],
                            np = 3)
    azimuth <- predict(fit, row, col) %>% round()

    i <- i[!u]
    zenith[!u] <- z[i][,]
    azimuth[!u] <- a[i][,]
  } else {
    zenith <- z[i][,]
    azimuth <- a[i][,]
  }

  data.frame(zenith, azimuth)
}
