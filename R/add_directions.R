.add_directions <- function(sampling_indices, z, a, in_radians = TRUE) {
  sampling_indices <- extract_dn(
    c(z, a),
    sampling_indices,
    use_window = FALSE
  )
  angular_columns <- ncol(sampling_indices) - c(1, 0)
  names(sampling_indices)[angular_columns] <- c("zen", "azi")
  if (in_radians) {
    sampling_indices$zen <- .degree2radian(sampling_indices$zen)
    sampling_indices$azi <- .degree2radian(sampling_indices$azi)
  }
  sampling_indices
}

#' Add angular directions to the data frame
#'
#' @description
#' Enrich a data frame containing cell indices (`"col"`, `"row"`) with zenith
#' and azimuth angles (`"zen"`, `"azi"`) corresponding to those cells,
#' preserving any additional columns.
#'
#' @details In \pkg{rcaiman}, a direction refers to a viewing direction defined
#'   by zenith and azimuth angles in radians. These angles describe the
#'   orientation of a unit vector ray originating at the camera.
#'
#'   Cell indices are passed to [extract_dn()] for direct lookup in `z` and `a`.
#'
#' @note This function assumes that `z`, `a`, and `sampling_indices` are tied
#'   to the same underlying raster geometry.
#'
#'   If any extracted zenith value is `NA`, the function stops with an error.
#'
#' @param z [terra::SpatRaster-class] generated with [zenith_image()].
#' @param a [terra::SpatRaster-class] generated with [azimuth_image()].
#' @param sampling_indices `data.frame` with columns `"col"` and `"row"`,
#'   representing cell indices numbered from 1. It may contain additional
#'   columns.
#'
#' @return `data.frame` equal to `sampling_indices` with added columns `"zen"`
#'   (zenith angle, radians) and `"azi"` (azimuth angle, radians).
#'
#'
#'
#' @export
#'
#' @examples
#' z <- zenith_image(10, equidistant_coef)
#' a <- azimuth_image(z)
#' sampling_indices <- data.frame(
#'   col = c(5, 6, 8),
#'   row = c(5, 5, 6)
#' )
#' add_directions(sampling_indices, z, a)
add_directions <- function(sampling_indices, z, a) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_columns(
    sampling_indices,
    required_columns = c("col", "row"),
    reserved_columns = c("zen", "azi")
  )

  sampling_indices <- .add_directions(
    sampling_indices,
    z,
    a,
    in_radians = TRUE
  )
  if (any(is.na(sampling_indices$zen)))
    stop("Some sampling points fall outside the valid zenith angle ",
         "domain (`zen` is `NA`).")
  sampling_indices
}
