#' Obtain row and col numbers from zenith and azimuth angles
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams zenith_image
#' @param zenith,azimuth Numeric vector. Zenith or azimuth angle in degrees.
#'
#'
#' @export
#'
#' @return An object of the class _data.frame_.
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#' row_col_from_zenith_azimuth(z, a, 45, 270)
row_col_from_zenith_azimuth <- function(z, a, zenith, azimuth) {
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(z, a)
  stopifnot(length(zenith) == length(azimuth))

  size <- 100
  sky_points <- expand.grid(row = seq(1, nrow(z), length.out = size) %>%
                              round(),
                            col = seq(1, nrow(z), length.out = size) %>%
                              round())
  sky_points <- extract_dn(z, sky_points, use_window = FALSE)
  sky_points <- sky_points[!is.na(sky_points[,3]), c("row", "col")]
  sky_points <- extract_rel_radiance(z, z, a, sky_points,
                                     use_window = FALSE)$sky_points

  fit <- spatial::surf.ls(x = sky_points[, "a"],
                          y = sky_points[, "z"],
                          z = sky_points[, "row"],
                          np = 6)
  row <- predict(fit, azimuth, zenith) %>% round()

  fit <- spatial::surf.ls(x = sky_points[, "a"],
                          y = sky_points[, "z"],
                          z = sky_points[, "col"],
                          np = 6)
  col <- predict(fit, azimuth, zenith) %>% round()

  data.frame(row, col)
}
