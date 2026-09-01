#' Reproject hemispherical image between radial projection models
#'
#' @description
#' Recompute raster values after changing the radial projection model that
#' defines pixel radial distances.
#'
#' @details
#' The function samples all the pixels of the cicular image with
#' [raster_grid_in_circle()], reprojects with [reproject_points()], and produce
#' and image with [interpolate_planar()].
#'
#' Arguments `z`, `a`, and `m` describe `r`, while `projection_coef` and
#' `diameter` describe the target, and `k`, `p`, and `search_radius` control
#' the interpolation.
#'
#' @param r numeric [terra::SpatRaster-class]
#'
#' @inheritParams compute_canopy_openness
#' @inheritParams add_directions
#' @inheritParams reproject_points
#' @inheritParams interpolate_planar
#'
#' @return numeric [terra::SpatRaster-class] with the same number of layers as
#'   `r`, reprojected under the target radial projection model.
#'
#' @export
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r1 <- reproject_image(r, z, a, m, equidistant_coef, ncol(r))
#' r2 <- reproject_image(r, z, a, m, orthographic_coef, ncol(r))
#'
#' plot(r)
#' plot(r1)
#' plot(r2)
#' }
reproject_image <- function(r, z, a, m,
                            projection_coef,
                            diameter,
                            k = 1,
                            p = 1,
                            search_radius = 10,
                            parallel = TRUE,
                            cores = NULL,
                            logical = TRUE,
                            leave_free = 0) {

  .assert_spatraster(r)
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_logical_mask(m)
  .assert_same_geom(z, r)
  .assert_same_geom(z, a)
  .assert_same_geom(m, r)

  # sample grid (angular domain)
  full_grid <- raster_grid_in_circle(z, a, 1, m)

  .reproject_layer <- function(layer) {

    # extract values
    pts <- extract_dn(layer, full_grid, use_window = FALSE)

    # reproject coordinates (new model)
    pts <- reproject_points(pts,
                            projection_coef = projection_coef,
                            diameter = diameter)

    # interpolate to raster
    interpolate_planar(pts, layer,
                       col_id = names(layer),
                       engine = "cpp",
                       search_radius = search_radius,
                       parallel = parallel,
                       cores = cores,
                       logical = logical,
                       leave_free = leave_free)
  }
  layer_names <- names(r)
  out <- Map(.reproject_layer, as.list(r))
  out <- terra::rast(out)
  names(out) <- layer_names
  out[!m] <- 0
  out
}
