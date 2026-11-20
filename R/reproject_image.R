#' Reproject an image to a radial destination
#'
#' @description
#' Resample image values onto a radial destination grid using planar
#' inverse-distance interpolation.
#'
#' @param r numeric [terra::SpatRaster-class] with one or more layers.
#'   Apply radiometric corrections to these original data before reprojection.
#' @param d numeric [terra::SpatRaster-class] describing the source directions,
#'   as returned by [direction_image()]. Must contain uniquely named `zenith`
#'   and `azimuth` layers in degrees and have the same grid as `r`.
#' @param m optional logical [terra::SpatRaster-class] with one layer, no `NA`,
#'   and the same grid as `r`. Only source cells where `m` is `TRUE` supply
#'   interpolation samples. If `NULL`, all source cells are considered.
#' @param dim numeric vector of length two with positive integer-like values.
#'   Destination dimensions in the order `c(width, height)`.
#'
#' @inheritParams reproject_points
#' @inheritParams interpolate_planar
#'
#' @details
#' By default, the destination is centered on zenith with North at the top.
#' Use `sensor_rotation`, `azimuth_rotation`, and `elevation` to orient the
#' destination as in [direction_image()]. Source orientation is already
#' described by `d`. Supply the same
#' `dim` and radial `geometry` that would construct its synthetic internal
#' geometry with [set_internal_geometry()]. The source can be rectangular,
#' cropped, or inclined, and its directions can come from either lens model.
#'
#' To align acquisitions on a common representation, adjust their source
#' orientations through `d` and keep the destination geometry, dimensions, and
#' orientation the same. See "Describing an imperfectly levelled acquisition"
#' in [direction_image()] for an example starting from an upward-looking
#' position.
#'
#' Source pixel centers are mapped to continuous destination coordinates with
#' [reproject_points()]. Their values are interpolated onto the destination
#' grid with [interpolate_planar()]. Each band uses its own available samples;
#' missing values and directions do not supply observations. The number and
#' names of bands are preserved.
#'
#' `search_radius` is measured in destination pixels. Together with `k`, `p`,
#' and `rigorous`, it controls interpolation. A source mask selects samples;
#' it does not become a destination mask. Interpolation can fill gaps and
#' extend values near the edges of the sampled region. Its reach should not
#' be interpreted as the observed coverage. Destination cells without enough
#' neighbors, or outside the radial angular limit, contain `NA`.
#'
#' Reprojection resamples image values and is generally not reversible.
#'
#' @return numeric [terra::SpatRaster-class] with dimensions `dim`, the same
#'   layers as `r`, and the standardized grid of [set_rcaiman_geometry()].
#'
#' @seealso [reproject_points()], [direction_image()]
#'
#' @export
#' @examples
#' source_geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(20, 15),
#'   parameters = 1 / (pi / 3),
#'   normalization_radius = 25,
#'   normalization_angle = 60
#' )
#' g <- set_internal_geometry(c(40, 30), source_geometry)
#' d <- direction_image(g, elevation = 32.5)
#'
#' # A synthetic radiometrically corrected band for this example.
#' r <- g[[1]] * 0 + 1
#' names(r) <- "radiance"
#' target_geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(30, 30),
#'   parameters = equidistant_coef,
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' standard <- reproject_image(
#'   r, d,
#'   dim = c(60, 60),
#'   geometry = target_geometry,
#'   search_radius = 1,
#'   parallel = FALSE
#' )
#' plot(standard)
#'
#' # A destination inclined toward south, with its sensor in portrait.
#' inclined <- reproject_image(
#'   r, d,
#'   dim = c(60, 60),
#'   geometry = target_geometry,
#'   sensor_rotation = 90,
#'   elevation = 32.5,
#'   search_radius = 1,
#'   parallel = FALSE
#' )
reproject_image <- function(r,
                            d,
                            dim,
                            geometry,
                            sensor_rotation = 0,
                            azimuth_rotation = 0,
                            elevation = 90,
                            m = NULL,
                            k = 1,
                            p = 1,
                            search_radius = 10,
                            rigorous = FALSE,
                            engine = "cpp",
                            parallel = TRUE,
                            cores = NULL,
                            logical = TRUE,
                            leave_free = 0) {

  .assert_spatraster(r)
  .assert_spatraster(d)
  .assert_same_geom(r, d)
  if (anyDuplicated(names(d)) ||
      !all(c("zenith", "azimuth") %in% names(d))) {
    stop("`d` must contain uniquely named `zenith` and `azimuth` layers.",
         call. = FALSE)
  }
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, r)
  }
  .check_radial_destination(geometry)
  target_g <- set_internal_geometry(dim, geometry)
  template <- target_g[[1]]
  outside <- is.na(terra::values(template, mat = FALSE))
  d <- d[[c("zenith", "azimuth")]]
  directions <- terra::values(d)
  selected <- !is.na(directions[, "zenith"]) &
    (!is.na(directions[, "azimuth"]) | directions[, "zenith"] %in% c(0, 180))
  if (!is.null(m)) selected <- selected & terra::values(m, mat = FALSE)
  cells <- which(selected)
  indices <- terra::rowColFromCell(d, cells)
  points <- add_directions(data.frame(col = indices[, 2], row = indices[, 1]), d = d)
  points <- reproject_points(points, geometry, sensor_rotation, azimuth_rotation, elevation)
  selected <- is.finite(points$x) & is.finite(points$y_down)
  cells <- cells[selected]
  points <- points[selected, c("x", "y_down"), drop = FALSE]

  layers <- lapply(seq_len(terra::nlyr(r)), function(i) {
    samples <- points
    samples$dn <- terra::values(r[[i]], mat = FALSE)[cells]
    out <- interpolate_planar(
      samples,
      r = template,
      k = k,
      p = p,
      search_radius = search_radius,
      rigorous = rigorous,
      engine = engine,
      parallel = parallel,
      cores = cores,
      logical = logical,
      leave_free = leave_free
    )
    values <- terra::values(out, mat = FALSE)
    values[outside] <- NA_real_
    terra::setValues(out, values)
  })
  out <- terra::rast(layers)
  names(out) <- names(r)
  out
}
