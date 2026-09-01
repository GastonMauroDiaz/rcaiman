#' Update raster coordinates after radial reprojection
#'
#' @description
#' Recompute raster coordinates of points after changing the radial projection
#' model that defines their radial distances.
#'
#' @details
#' \pkg{rcaiman} represents hemispherical images using a radial projection
#' system, where each point is defined by its angular coordinates (zenith and
#' azimuth).
#'
#' The radial projection model defines how zenith angles are transformed into
#' normalized radial distances (see [calc_relative_radius()]). These distances
#' are then scaled to a planar grid and discretized into pixel indices according
#' to `diameter`.
#'
#' Reprojection consists of modifying the radial distance of each point while
#' preserving its angular coordinates.
#'
#' @param sampling_points data.frame with at least columns `"row"` and `"col"`.
#'   If `use_existing_angles = TRUE`, it must also contain `"zen"` and `"azi"`
#'   (see [add_directions()]).
#'
#' @inheritParams calc_relative_radius
#' @inheritParams zenith_image
#'
#' @return data.frame with updated `row` and `col` corresponding to the
#'   reprojected coordinates.
#'
#' @export
#' @examples
#' z <- zenith_image(1000, equal_area_coef)
#' a <- azimuth_image(z)
#' pts <- fibonacci_points(z, a, 10)
#' pts2 <- reproject_points(pts, equidistant_coef, ncol(z))
#'
#' plot(!is.na(z))
#' points(pts$col, nrow(z) - pts$row, col = "blue", pch = 10)
#' points(pts2$col, nrow(z) - pts2$row, col = "darkgreen", pch = 10)
reproject_points <- function(sampling_points,
                             projection_coef,
                             diameter) {

  .check_columns(sampling_points,
                 required_columns = c("zen", "azi"))
  .check_vector(diameter, "even_integerish", 1, sign = "positive")

  zen <- sampling_points$zen
  azi <- sampling_points$azi

  R <- .calc_relative_radius(zen, projection_coef)

  theta <- azi + pi/2
  r <- R * (pi/2)

  x <- r * cos(theta)
  y <- r * sin(theta)

  n <- diameter
  res <- pi / n

  # col <- floor((x + pi/2) / res) + 1
  # row <- floor((pi/2 - y) / res) + 1
  col <- (x + pi/2) / res + 1
  row <- (pi/2 - y) / res + 1

  sampling_points$row <- row
  sampling_points$col <- col

  sampling_points
}
