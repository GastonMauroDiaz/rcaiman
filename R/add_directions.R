.add_directions <- function(sampling_points, z, a, in_radians = TRUE) {
  sampling_points <- extract_dn(c(z, a), sampling_points, use_window = FALSE)
  sp_names <- gsub("Zenith image", "zen", names(sampling_points))
  sp_names <- gsub("Azimuth image", "azi", sp_names)
  names(sampling_points) <- sp_names
  if (in_radians) {
    sampling_points$zen <- .degree2radian(sampling_points$zen)
    sampling_points$azi <- .degree2radian(sampling_points$azi)
  }
  sampling_points
}

#' Add angular directions to the data frame
#'
#' @description
#' Enrich a data frame containing continuous image coordinates (see
#' [xy_from_row_col()]) with zenith and azimuth angles (`"zen"`, `"azi"`)
#' corresponding to those locations, preserving any additional columns.
#'
#' @details In \pkg{rcaiman}, a direction refers to a viewing direction defined
#'   by zenith and azimuth angles in radians. These angles describe the
#'   orientation of a unit vector ray originating at the camera.
#'
#' @note This function assumes that `z`, `a`, and `sampling_points` are tied to
#'   the same underlying raster geometry.
#'
#'   If any extracted zenith value is `NA`, the function stops with an error.
#'
#' @param z [terra::SpatRaster-class] generated with [zenith_image()].
#' @param a [terra::SpatRaster-class] generated with [azimuth_image()].
#' @param sampling_points `data.frame` with columns `row` and `col` (raster
#'   coordinates). It may contain additional columns.
#'
#' @return `data.frame` equal to `sampling_points` with added columns `"zen"`
#'   (zenith angle, radians) and `"azi"` (azimuth angle, radians).
#'
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#'
#' z <- zenith_image(spec$horizon_radius * 2, spec$parameters)
#' a <- azimuth_image(z)
#'
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' head(add_directions(sky_points, z, a))
#' }
add_directions <- function(sampling_points, z, a) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_columns(
    sampling_points,
    required_columns = c("row", "col"),
    reserved_columns = c("zen", "azi")
  )

  sampling_points <- .add_directions(sampling_points, z, a, in_radians = TRUE)
  if(any(is.na(sampling_points$z)))
    stop("Some sampling points fall outside the valid zenith angle ",
         "domain (`zen` is `NA`).")
  sampling_points
}
