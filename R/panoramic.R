#' Transform hemispherical data through a panoramic domain
#'
#' @description
#' Reparameterize hemispherical angular coordinates into a panoramic
#' working domain, operate in this domain following a planar
#' logic, and map the resulting values back into hemispherical image domain.
#'
#' @details
#' These functions transform hemispherical angular coordinates into a temporary
#' panoramic domain.
#'
#' Continuous image coordinates are derived directly from zenith and azimuth
#' angles, without selecting cells. After padding, `x` increases with azimuth
#' and `y_down` increases with zenith angle. One coordinate unit corresponds
#' to `angle_width` degrees. The unpadded origin is `zen = 0`, `azi = 0`.
#'
#' This introduces artificial discontinuities that are not
#' present in hemispherical geometry:
#'
#' * azimuth becomes discontinuous at the 0 deg / 360 deg seam,
#' * and directions near the zenith become artificially separated in the
#'   rectangular domain despite converging geometrically toward a single point.
#'
#' To reduce these effects during computation, the panoramic domain is extended
#' with auxiliary padding regions:
#'
#' * lateral padding in both sides of the azimuth seam,
#' * and reflected padding across the zenith boundary.
#'
#' @section Functions:
#' \describe{
#'   \item{`to_panoramic()`}{
#'   Transform hemispherical angular coordinates into a temporary panoramic
#'   raster domain with or without padding regions.
#'   }
#'
#'   \item{`rem_panoramic_padding()`}{
#'   Remove the padding regions from a panoramic raster.
#'   }
#'
#'   \item{`from_panoramic()`}{
#'   Map values from panoramic geometry back into hemispherical image geometry
#'   using zenith and azimuth images.
#'   }
#' }
#'
#' @param sampling_points `data.frame` containing numeric columns `"zen"` and
#'   `"azi"` in radians, as in the output of [add_directions()] with
#'   `in_radians = TRUE`, or [add_directions_from_model()]. Existing `"x"` and
#'   `"y_down"` columns are replaced by panoramic coordinates. Derived columns
#'   tied to the source image frame are removed with [keep_coord()]. Other
#'   columns are preserved, including the source directions.
#' @param angle_width numeric vector of length one. Angular resolution in
#'   degrees.
#' @param padding_deg integer-like numeric vector of length one. Nonnegative
#'   width of the padding region in degrees. Use zero for no padding.
#' @param pano single-layer [terra::SpatRaster-class] in panoramic geometry.
#'   Remove padding before passing this raster to `from_panoramic()`.
#' @param z [terra::SpatRaster-class] generated with [zenith_image()].
#' @param a [terra::SpatRaster-class] generated with [azimuth_image()]
#'
#' @return
#' `to_panoramic()` returns a list with named elements:
#' \describe{
#'   \item{`sampling_coord`}{A `data.frame` with continuous `x`, `y_down`
#'     coordinates and preserved frame-independent columns. Original rows
#'     appear first, followed by auxiliary top, left, and right padding rows.}
#'   \item{`raster`}{An empty, single-layer [terra::SpatRaster-class] with
#'     geometry set by [set_rcaiman_geometry()], for planar interpolation.}
#' }
#' `rem_panoramic_padding()` returns a raster without the padding regions.
#' `from_panoramic()` returns a raster with the geometry of `z`, populated
#' from the panoramic cells containing the corresponding directions. At the
#' outer boundaries, indices are clamped to the nearest valid cell. Cells
#' without a zenith value remain `NA`.
#'
#' @seealso [fisheye_to_pano()]
#'
#' @name panoramic
#' @rdname panoramic
#' @examples
#' z <- zenith_image(30, equidistant_coef)
#' a <- azimuth_image(z)
#' sampling_indices <- raster_grid_in_circle(z, a, spacing = 3)
#' sampling_indices$dn <- cos(sampling_indices$zen)
#'
#' angle_width <- 5
#' padding_deg <- 10
#' pano <- to_panoramic(sampling_indices, angle_width, padding_deg)
#' interpolated <- interpolate_planar(
#'   pano$sampling_coord, r = pano$raster, search_radius = 20,
#'   parallel = FALSE
#' )
#' unpadded <- rem_panoramic_padding(interpolated, angle_width, padding_deg)
#' hemi <- from_panoramic(unpadded, z, a, angle_width)
#' plot(hemi)
to_panoramic <- function(sampling_points,
                         angle_width = 1,
                         padding_deg = 15) {

  .check_columns(
    sampling_points,
    required_columns = c("zen", "azi"),
    reserved_columns = NULL
  )

  .check_vector(angle_width, "numeric", 1, sign = "positive")
  .check_vector(padding_deg, "integerish", 1, sign = "nonnegative")

  sampling_coord <- sampling_points
  # Work in degrees for padding, then scale to the raster's unit geometry.
  sampling_coord$y_down <- sampling_coord$zen * 180 / pi + padding_deg
  sampling_coord$x <- sampling_coord$azi * 180 / pi + padding_deg
  sampling_coord <- keep_coord(sampling_coord)

  # ---------------------------------------------------------------------------
  # Padding
  # ---------------------------------------------------------------------------

  top <- sampling_coord[sampling_coord$zen * 180 / pi < padding_deg, ]
  top$y_down <- padding_deg - (top$y_down - padding_deg)
  sampling_coord <- rbind(sampling_coord, top)

  left <- sampling_coord[sampling_coord$azi * 180 / pi < padding_deg, ]
  left$x <- left$x + 360

  right <- sampling_coord[sampling_coord$azi * 180 / pi > (360 - padding_deg), ]
  right$x <- right$x - 360

  sampling_coord <- rbind(sampling_coord, left, right)

  # ---------------------------------------------------------------------------
  # Raster template
  # ---------------------------------------------------------------------------

  ncol <- (360 + 2 * padding_deg) / angle_width
  nrow <- (90 + padding_deg ) / angle_width

  msn <- "This specific combination of `angle_width` and `padding_deg` cannot be used"
  res <- try(.check_vector(ncol, "integerish"))
  if(inherits(res, "try-error")) {
    stop(msn)
  }
  res <- try(.check_vector(nrow, "integerish"))
  if(inherits(res, "try-error")) {
    stop(msn)
  }

  r <- terra::rast(
    nrows = nrow,
    ncols = ncol,
    xmin = -padding_deg,
    xmax = 360 + padding_deg,
    ymin = 0,
    ymax = 90 + padding_deg
  )

  r <- set_rcaiman_geometry(r)

  sampling_coord$y_down <- sampling_coord$y_down / angle_width
  sampling_coord$x <- sampling_coord$x / angle_width

  list(
    sampling_coord = sampling_coord,
    raster = r
  )
}

#' @rdname panoramic
rem_panoramic_padding <- function(pano,
                                  angle_width = 1,
                                  padding_deg = 15) {
  .assert_spatraster(pano)
  .check_vector(angle_width, "numeric", 1, sign = "positive")
  .check_vector(padding_deg, "integerish", 1, sign = "nonnegative")

  mar <- padding_deg / angle_width
  xmn <- mar
  xmx <- mar + 360 / angle_width
  ymx <- 90 / angle_width
  ymn <- 0
  e <- terra::ext(xmn, xmx, ymn, ymx)
  pano <- terra::crop(pano, e)
  set_rcaiman_geometry(pano)
}

#' @rdname panoramic
from_panoramic <- function(pano,
                           z,
                           a,
                           angle_width = 1) {

  .assert_single_layer(pano)
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(angle_width, "numeric", 1, sign = "positive")

  idx <- which(!is.na(terra::values(z)))

  zen_deg <- terra::values(z)[idx]
  azi_deg <- terra::values(a)[idx]

  # Correction (2026-09-04): these ratios are continuous image coordinates,
  # not one-based cell indices. Select containing cells without rounding them
  # as indices; retain the existing clamping at the outer boundaries.
  pano_indices <- indices_from_coord(
    x = azi_deg / angle_width, y_down = zen_deg / angle_width
  )

  pano_indices$col[pano_indices$col < 1] <- 1
  pano_indices$col[pano_indices$col > ncol(pano)] <- ncol(pano)

  pano_indices$row[pano_indices$row < 1] <- 1
  pano_indices$row[pano_indices$row > nrow(pano)] <- nrow(pano)

  pano_cells <- terra::cellFromRowCol(
    pano,
    row = pano_indices$row,
    col = pano_indices$col
  )

  vals <- terra::values(pano)[pano_cells]

  out <- z
  terra::values(out) <- NA_real_
  terra::values(out)[idx] <- vals

  out
}
