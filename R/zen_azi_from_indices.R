#' Map between zenith-azimuth angles and cell indices
#'
#' @description
#' Convert between angular coordinates on a hemispherical image (`"zen"`, `"azi"`)
#' and cell indices (`"col"`, `"row"`). See
#' [y_down_from_y()] for coordinate conventions.
#'
#' @section Functions:
#' \describe{
#'   \item{`indices_from_zen_azi()`}{Return a `data.frame` with columns
#'     `"col"` and `"row"`. Each row corresponds to the input direction at the
#'     same position. Duplicated cells and input order are preserved.}
#'   \item{`zen_azi_from_indices()`}{Return a `data.frame` with columns `"zen"` and
#'     `azi` (radians) corresponding to the input cell indices.}
#' }
#'
#' @details
#' **`(zen, azi)` from `indices`.**
#' Angles are  obtained by direct lookup on `z` and `a` using the corresponding
#' raster cells. If the lookup fails (e.g., because the queried cells fall
#' outside the calibrated lens footprint), temporal angular images that
#' extrapolate beyond the horizon are reconstructed using the radial projection
#' model stored in `z` and the orientation attribute stored in `a`. Values are
#' translated from deg to radians.
#'
#' **`indices` from `(zen, azi)`.**
#' Zenith angles are first converted to relative radius using the radial
#' projection model stored in `z`, and azimuth angles are adjusted using the
#' orientation stored in `a`. The resulting polar coordinates are then
#' transformed to x-y coordinates and assigned directly to cells. Directions
#' outside the raster extent return `NA` indices. Directions falling in the
#' same cell remain as separate rows, and input order is preserved.
#'
#' @inheritParams add_directions
#' @inheritParams zenith_image
#'
#' @param zen,azi numeric vectors of equal length. Angles in radians.
#' @param col,row numeric vectors of equal length. Must be positive and
#'   interpreted as integer. Cell indices. See [y_down_from_y()].
#'
#' @return See *Functions*
#'
#' @name zen_azi_from_indices
#' @rdname zen_azi_from_indices
#' @aliases indices_from_zen_azi
#'
#' @export
#'
#' @examples
#' z <- zenith_image(1000, equidistant_coef)
#' a <- azimuth_image(z)
#'
#' indices <- indices_from_zen_azi(
#'   z,
#'   a,
#'   zen = c(30, 60) * pi / 180,
#'   azi = c(90, 270) * pi / 180
#' )
#' indices
#'
#' ang <- zen_azi_from_indices(
#'   z,
#'   a,
#'   col = indices$col,
#'   row = indices$row
#' )
#' ang
zen_azi_from_indices <- function(z, a, col, row) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)

  if (length(col) != length(row))
    stop("`col` and `row` must have same length.")

  .check_vector(col, "integerish", sign = "positive")
  .check_vector(row, "integerish", sign = "positive")

  if (max(col) > ncol(z))
    stop(sprintf("`col` must be up to %s", ncol(z)))

  if (max(row) > nrow(z))
    stop(sprintf("`row` must be up to %s", nrow(z)))

  sampling_indices <- data.frame(col, row)
  ang <- try(
    .add_directions(
      sampling_indices,
      z,
      a,
      in_radians = TRUE
    ),
    silent = TRUE
  )

  if (inherits(ang, "try-error")) {
    projection_coef <- attr(z, "projection_coef")
    diameter <- ncol(z)

    z <- .zenith_image(
      diameter,
      projection_coef,
      na_below_horizon = FALSE
    )
    a <- .azimuth_image(z, attr(a, "orientation"))

    ang <- .add_directions(
      sampling_indices,
      z,
      a,
      in_radians = TRUE
    )
  }

  ang[, c("zen", "azi")]
}


#' @rdname zen_azi_from_indices
#' @export
indices_from_zen_azi <- function(z, a, zen, azi) {
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

  # Project to relative radius using projection model
  R <- calc_relative_radius(
    zen,
    attr(z, "projection_coef"),
    units = "radian"
  )

  # Convert to projected Cartesian coordinates used by the image
  orientation <- attr(a, "orientation") * pi / 180
  pol <- data.frame( # col names required by geometry package.
    theta = azi + orientation + pi / 2, # This is not zenith angle, it is azimuth angle
    z = R * pi / 2 # This is not zenith angle, this is distance to center for a circle with radius = pi / 2
  )
  cart <- geometry::pol2cart(as.matrix(pol))

  # Use the angular raster geometry to locate each point in cell space.
  ztmp <- terra::deepcopy(z)
  terra::ext(ztmp) <- terra::ext(
    -pi / 2,
    pi / 2,
    -pi / 2,
    pi / 2
  )

  cell <- terra::cellFromXY(ztmp, cart)
  indices <- terra::rowColFromCell(ztmp, cell)
  indices <- as.data.frame(indices)
  names(indices) <- c("row", "col")
  indices[, c("col", "row")]
}
