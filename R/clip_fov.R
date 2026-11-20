#' Restrict an internal geometry to a field of view
#'
#' @description
#' Restrict the internal geometry of an imaging system to a supplied field of
#' view. For example, represent the field of view of a narrower-angle system on
#' the image grid of a wider-angle system. Directions outside the supplied
#' field of view are set to `NA`. The grid and its dimensions are preserved.
#'
#' @details
#' Each pixel-center direction is tested on the sphere, including the contour
#' boundary, with the numerical tolerance of [select_directions_in_fov()].
#' See [intersect_fov()] for the interpretation of contour
#' order and edges. The contour and the internal geometry must use the same
#' internal reference frame. No angles or image values are interpolated.
#'
#' A pixel with a missing angle is excluded from both layers. Existing missing
#' regions are preserved. To obtain a logical mask from the returned geometry,
#' use `!is.na(g$off_axis)`.
#'
#' @param g numeric [terra::SpatRaster-class] with layers
#'   `off_axis` and `internal_azimuth`, in that order, in degrees, as returned
#'   by [set_internal_geometry()].
#' @inheritParams select_directions_in_fov
#'
#' @return numeric [terra::SpatRaster-class] with the grid and layer names of
#'   `g`. Both angles are unchanged inside the field of view
#'   and `NA` elsewhere. The input geometry is not modified.
#'
#' @seealso [extract_fov()], [intersect_fov()],
#'   [select_directions_in_fov()], [rem_unselected()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(20, 15),
#'   normalization_radius = 18,
#'   normalization_angle = 90
#' )
#' wide_g <- set_internal_geometry(c(40, 30), geometry)
#' geometry$normalization_radius <- 30
#' narrow_g <- set_internal_geometry(c(40, 30), geometry)
#' narrow_fov <- extract_fov(narrow_g)
#' restricted_wide_g <- clip_fov(wide_g, narrow_fov)
#' plot(restricted_wide_g)
#' m <- !is.na(restricted_wide_g$off_axis)
clip_fov <- function(g, fov) {
  angles <- .internal_coverage_angles(g, "g")
  polygon <- .spherical_field_of_view(fov, "fov")
  directions <- as.data.frame(angles * pi / 180)
  selected <- .directions_in_field_of_view(directions, polygon)
  angles[!selected, ] <- NA_real_
  result <- terra::rast(g)
  terra::values(result) <- angles
  result
}
