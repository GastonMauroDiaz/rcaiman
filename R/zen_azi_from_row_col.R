#' Map between zenith–azimuth angles and raster coordinates
#'
#' @description
#' Convert between angular coordinates on a hemispherical image (`zen`, `azi`)
#' and continuous image coordinates (see [xy_from_col_row()]).
#'
#' @section Functions:
#' \describe{
#'   \item{`row_col_from_zen_azi`}{Return a `data.frame` with columns `col`
#'   and `row` corresponding to the input angular coordinates.}
#'   \item{`zen_azi_from_row_col`}{Return a `data.frame` with columns `zen`
#'   and `azi` (radians) corresponding to the input raster coordinates.}
#' }
#'
#' @details
#' **`(zen, azi)` from `(row, col)`.**
#' Angles are obtained by direct lookup on `z` and `a` using
#' [terra::cellFromRowCol()]. If the lookup fails (e.g., because the queried
#' cells fall outside the calibrated lens footprint), synthetic angular images
#' are reconstructed using the radial projection model stored in `z`
#' and the orientation attribute stored in `a`. Values are
#' translated from deg to radians.
#'
#' **`(row, col)` from `(zen, azi)`.**
#' Zenith angles are first converted to relative radius using the radial
#' projection model stored in `z`. The resulting polar coordinates are then
#' transformed to Cartesian coordinates and rasterized to identify the
#' corresponding pixel indices.
#'
#' If this direct mapping fails, a fallback approximation is used: a sparse set
#' of sky points generated with [fibonacci_points()] is enriched with angular
#' coordinates and used to fit two local least-squares surfaces
#' ([spatial::surf.ls()]) predicting `row` and `col` from `azi` and `zen`.
#'
#' @inheritParams add_directions
#' @inheritParams zenith_image
#'
#' @param zen,azi numeric vectors of equal length. Angles in radians.
#' @param row,col numeric vectors of equal length. Continuous image coordinates
#'   (see [xy_from_col_row()]).
#'
#' @return See *Functions*
#'
#' @name zen_azi_from_row_col
#' @rdname zen_azi_from_row_col
#' @aliases row_col_from_zen_azi
#'
#' @export
#'
#' @examples
#' z <- zenith_image(1000, equidistant_coef)
#' a <- azimuth_image(z)
#'
#' rc <- row_col_from_zen_azi(z, a, zen = c(30, 60), azi = c(90, 270))
#' rc
#'
#' ang <- zen_azi_from_row_col(z, a, row = rc$row, col = rc$col)
#' ang
zen_azi_from_row_col <- function(z, a, row, col) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  if (length(row) != length(col))
    stop("`row` and `col` must have same length.")
  .check_vector(row, "numeric", sign = "positive")
  .check_vector(col, "numeric", sign = "positive")
  if (max(row) > ncol(z))
    stop(sprintf("`row` cannot be greater than %s", ncol(z)))
  if (max(col) > ncol(z))
    stop(sprintf("`col` cannot be greater than %s", ncol(z)))
  if (length(row) != length(col))
    stop("`row` and `col` must have same length.")

  ang <- try(.add_directions(data.frame(row, col), z, a, in_radians = TRUE))
  if(inherits(ang, "try-error")) {
    projection_coef <- attr(z, "projection_coef")
    diameter <- ncol(z)
    z <- .zenith_image(diameter, projection_coef, na_below_horizon = FALSE)
    a <- .azimuth_image(z, attr(a, "orientation")) #it takes NA from z
    ang <- .add_directions(data.frame(row, col), z, a, in_radians = TRUE)
  }
  ang[, c("zen", "azi")]
}

#' @rdname zen_azi_from_row_col
#' @export
row_col_from_zen_azi <- function(z, a, zen, azi) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(zen, "numeric", sign = "nonnegative")
  .check_vector(azi, "numeric", sign = "nonnegative")
  if (length(zen) != length(azi))
    stop("`zen` and `azi` must have same length.")
  if (any(zen > pi))
    stop("`zen` cannot be greater than pi radians (180 deg).")
  if (any(azi > 2 * pi))
    stop("`azi` cannot be greater than 2pi radians (360 deg).")

  .row_col_from_zen_azi <- function(z, a, zen, azi) {
    # project to relative radius using projection model
    R <- calc_relative_radius(zen, attr(z, "projection_coef"), units = "radian")

    # convert to projected Cartesian coordinates used by the image
    pol <- data.frame(
      theta = azi + pi/2,
      z = R * pi/2
    )
    cart <- geometry::pol2cart(as.matrix(pol))
    p <- terra::vect(matrix(cart, ncol = 2))
    terra::crs(p) <- terra::crs(z)

    # rasterize to locate each point in pixel space
    ztmp <- terra::deepcopy(z)
    terra::ext(ztmp) <- terra::ext(-pi/2, pi/2, -pi/2, pi/2)

    r <- terra::rasterize(p, ztmp)
    r <- !is.na(r)
    idx <- cells(r)[r[]]
    row_col <- terra::rowColFromCell(ztmp, idx)
    row_col <- as.data.frame(row_col)
    names(row_col) <- c("row_index", "col_index")

    col_row_from_cell_indices(row_col$col_index, row_col$row_index)
  }

  row_col <- try(.row_col_from_zen_azi(z, a, zen, azi))
  if(inherits(row_col, "try-error")) {
    sampling_points  <- fibonacci_points(z, a, 10)

    fit <- spatial::surf.ls(x = sky_points[, "azi"],
                            y = sky_points[, "zen"],
                            z = sky_points[, "row"],
                            np = 6)
    row <- predict(fit, azi, zen)

    fit <- spatial::surf.ls(x = sky_points[, "azi"],
                            y = sky_points[, "zen"],
                            z = sky_points[, "col"],
                            np = 6)
    col <- predict(fit, azi, zen)
    row_col <- data.frame(row, col)
  }
  row_col
}
