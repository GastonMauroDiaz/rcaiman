#' Set rcaiman internal raster geometry
#'
#' @description
#' Apply the internal raster geometry convention used throughout
#' *rcaiman*: rasters are treated as image grids, not as spatial objects.
#'
#' @details
#' This function enforces a unit-resolution pixel grid with a
#' standardized extent `ext(0, ncol, 0, nrow)` and assigns a non-geographic CRS
#' as an explicit marker that no spatial semantics are intended. The `EPSG:7589`
#' was selected because it does not have datum, projection, or physical units.
#' This ensure consistent behavior across `terra`-based workflows. Using `crs =
#' NA` is intentionally avoided since many `terra` methods and external tools
#' assume a defined CRS and may otherwise fail, emit warnings, or apply implicit
#' defaults.
#'
#' @param r numeric [terra::SpatRaster-class].
#'
#' @return The same [terra::SpatRaster-class] with standardized extent and CRS.
#'
#' @export
set_rcaiman_geometry <- function(r) {
  .assert_spatraster(r)
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  terra::crs(r) <- "epsg:7589"
  r
}
