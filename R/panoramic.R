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
#' Coordinates are derived directly from zenith angle and azimuth angle.
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
#'   \item{`remove_panoramic_padding()`}{
#'   Remove the padding regions from a panoramic raster.
#'   }
#'
#'   \item{`from_panoramic()`}{
#'   Map values from panoramic geometry back into hemispherical image geometry
#'   using zenith and azimuth images.
#'   }
#' }
#'
#' @param sampling_points `data.frame` containing at least `"zen"` and `"azi"`,
#'   as in the output of [add_directions()]. Additional columns are preserved.
#' @param angle_width numeric vector of length one. Angular resolution in
#'   degrees.
#' @param padding_deg numeric vector of length one. Must be positive and
#'   interpreted as integer.  Width of the padding region in degrees.
#' @param pano single-layer [terra::SpatRaster-class] in panoramic geometry.
#' @param z [terra::SpatRaster-class] generated with [zenith_image()].
#' @param a [terra::SpatRaster-class] generated with [azimuth_image()]
#'
#' @return See *Functions*.
#'
#' @seealso [fisheye_to_pano()]
#'
#' @name panoramic
#' @rdname panoramic
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' rspec <- get_radiometry_spec(spec, "simple_method")
#'
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#'
#' # Manipulate image
#' sampling_points <- raster_grid_in_circle(z, a, 1)
#' sampling_points <- extract_dn(r, sampling_points, use_window = FALSE)
#'
#' angle_width <- 1
#' padding_deg <- 10
#' pano <- to_panoramic(sampling_points,
#'                      angle_width = angle_width,
#'                      padding_deg = padding_deg)
#' pano <- interpolate_planar(pano$sampling_points,
#'                            "Cyan",
#'                            pano$raster,
#'                            search_radius = 5)
#' plot(pano)
#' pano <- rem_panoramic_padding(pano, angle_width, padding_deg)
#' plot(pano)
#'
#' hemi <- from_panoramic(pano, z, a, angle_width = angle_width)
#' plot(hemi)
#' plot(r - hemi)
#'
#' # Work in panoramic and return to hemispherical
#'
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' sky_points <- extract_dn(r, sky_points, use_window = FALSE)
#' sky_points <- add_directions(sky_points, z, a)
#'
#' angle_width <- 1
#' padding_deg <- 40
#' pano <- to_panoramic(sky_points,
#'                      angle_width = angle_width,
#'                      padding_deg = padding_deg)
#'
#' ## IDW
#' idw <- interpolate_planar(pano$sampling_points,
#'                           "Cyan",
#'                           pano$raster,
#'                           search_radius = 20)
#' idw <- rem_panoramic_padding(idw, angle_width, padding_deg)
#' idw <- from_panoramic(idw, z, a, angle_width = angle_width)
#' plot(idw)
#' ## Triangulation
#' tri <- triangulate_planar(pano$sampling_points,
#'                           "Cyan",
#'                           pano$raster)
#' tri <- rem_panoramic_padding(tri, angle_width, padding_deg)
#' tri <- from_panoramic(tri, z, a, angle_width = angle_width)
#' plot(tri)
#'
#' ## Fit trend surface
#' ts <- fit_trend_surface(pano$sampling_points,
#'                         "Cyan",
#'                         pano$raster,
#'                         extrapolate = TRUE)
#' ts <- rem_panoramic_padding(ts, angle_width, padding_deg)
#' ts <- from_panoramic(ts, z, a, angle_width = angle_width)
#' plot(ts)
#' }
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

  sampling_points$row <- sampling_points$zen * 180 / pi + padding_deg
  sampling_points$col <- sampling_points$azi * 180 / pi + padding_deg

  # ---------------------------------------------------------------------------
  # Padding
  # ---------------------------------------------------------------------------

  top <- sampling_points[sampling_points$zen * 180 / pi < padding_deg, ]
  top$row <- padding_deg - (top$row - padding_deg)
  sampling_points <- rbind(sampling_points, top)

  left <- sampling_points[sampling_points$azi * 180 / pi < padding_deg, ]
  left$col <- left$col + 360

  right <- sampling_points[sampling_points$azi * 180 / pi > (360 - padding_deg), ]
  right$col <- right$col - 360

  sampling_points <- rbind(sampling_points, left, right)

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

  sampling_points$row <- sampling_points$row / angle_width
  sampling_points$col <- sampling_points$col / angle_width

  list(
    sampling_points = sampling_points,
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

  idx <- which(!is.na(terra::values(z)))

  zen_deg <- terra::values(z)[idx]
  azi_deg <- terra::values(a)[idx]

  col_pano <- round(azi_deg  / angle_width)
  row_pano <- round(zen_deg / angle_width)

  col_pano[col_pano < 1] <- 1
  col_pano[col_pano > ncol(pano)] <- ncol(pano)

  row_pano[row_pano < 1] <- 1
  row_pano[row_pano > nrow(pano)] <- nrow(pano)

  pano_cells <- terra::cellFromRowCol(
    pano,
    row = row_pano,
    col = col_pano
  )

  vals <- terra::values(pano)[pano_cells]

  out <- z
  terra::values(out) <- NA_real_
  terra::values(out)[idx] <- vals

  out
}
