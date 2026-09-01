#' Estimate sun angular coordinates
#'
#' Estimates the sun’s zenith and azimuth angles (deg) from a canopy
#' hemispherical photograph, using either direct detection of the solar disk or
#' indirect cues from the circumsolar region.
#'
#' This function can operate under two alternative assumptions for
#' estimating the sun position:
#'
#' \describe{
#'   \item{Veiled sun}{The solar disk is visible or partially obscured; the
#'   position is inferred from localized brightness peaks. It is also useful
#'   when a sunscreen is used, as in \insertCite{Kuusk2007;textual}{rcaiman}.}
#'   \item{Obscured sun}{The solar disk is not visible; the position is inferred
#'   from radiometric and spatial cues aggregated over the circumsolar region.}
#' }
#'
#' When `mode = "veiled"`, `seg` and `circumsolar_angle` are ignored.
#'
#' Estimates refer to positions above the horizon; therefore, estimated angles
#' may require further manipulation if the photograph was acquired under
#' crepuscular light.
#'
#' @note A scientific article presenting and validating this method is currently
#' under preparation.
#'
#' @param mode character vector of length one. Estimation mode:
#'   `"obscured"` (default) or `"veiled"`. See *Details*.
#' @param circumsolar_angle numeric vector of length one. Maximum angular
#'   radius (in degrees) used to define the circumsolar region.
#' @param seg Segmentation map of r, typically created with functions such as
#'   [equal_area_segmentation()], [sky_grid_segmentation()], [ring_segmentation()]
#'   or [sector_segmentation()], but any raster with integer segment labels is
#'   accepted.
#'
#' @inheritParams sample_sky_points
#' @inheritParams add_directions
#' @inheritParams compute_mean
#'
#' @seealso [optimize_sun_angles()]
#'
#' @returns Named numeric vector of length two, with names `z` and `a`,
#'   representing the sun’s zenith and azimuth angles (in degrees).
#'
#' @references \insertAllCited{}
#'
#' @export
#'
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
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' thr <- apply_by_direction(r, z, a, m)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' seg <- sky_grid_segmentation(z, a, 10)
#' sun_angles <- estimate_sun_angles(r, z, a, bin, seg,
#'                                   circumsolar_angle = 30,
#'                                   mode = "obscured")
#' row_col <- row_col_from_zen_azi(z, a,
#'                                 sun_angles["zen_deg"] * pi/180,
#'                                 sun_angles["azi_deg"] * pi/180)
#' plot(caim$Cyan)
#' points(row_col[1,2], nrow(caim) - row_col[1,1], col = "yellow",
#'        pch = 8, cex = 3)
#'
#'
#' sun_angles <- estimate_sun_angles(r, z, a, bin,
#'                                   mode = "veiled")
#' row_col <- row_col_from_zen_azi(z, a,
#'                                 sun_angles["zen_deg"] * pi/180,
#'                                 sun_angles["azi_deg"] * pi/180)
#' plot(caim$Cyan)
#' points(row_col[1,2], nrow(caim) - row_col[1,1], col = "yellow",
#'        pch = 8, cex = 3)
#' }
estimate_sun_angles <- function(r, z, a, bin, seg,
                                circumsolar_angle = 30,
                                mode = "obscured") {

  .this_requires_EBImage()

  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .assert_single_layer(r, allow_na = FALSE)
  .assert_same_geom(r, z)
  .assert_logical_mask(bin)
  .assert_same_geom(r, bin)
  .assert_choice(mode, c("obscured", "veiled"))

  if (mode == "obscured") {
    .assert_single_layer(seg)
    .assert_same_geom(r, seg)
    .check_vector(circumsolar_angle, "numeric", 1, sign = "positive")

    if (!is.numeric(circumsolar_angle) || length(circumsolar_angle) != 1) {
      stop("`circumsolar_angle` must be a numeric vector of length one.")
    }

    # Select cells with at least one extremely bright sky pixel
    seg[!bin] <- NA
    r <- compute_mean(r, seg, function(x) quantile(x, 0.99, na.rm = TRUE))
    m <- r >= quantile(unique(terra::values(r)), 0.9, na.rm = TRUE)
    m[is.na(m)] <- 0

    # Merge adjacent white segments
    labeled_m <- EBImage::bwlabel(as.array(m)[,,1])
    labeled_m <- terra::rast(labeled_m)
    terra::ext(labeled_m) <- terra::ext(r)
    terra::crs(labeled_m) <- terra::crs(r)
    # Calc membership to class "sun seed"
    .fun <- function(x) {
      x <- unique(x) # to count the cells instead of the pixels
      length(x)
    }
    size <- compute_mean(seg, labeled_m, .fun, return = "vector") %>%
      normalize_min_max()
    dn <- compute_mean(r, labeled_m, median, return = "vector") %>%
      normalize_min_max()
    if (any(is.nan(dn))) dn[] <- 1 #preventative programming
    if (any(is.nan(size))) size[] <- 1 #preventative programming
    membership_posibility <- size*0.25 + dn*0.75

    # Find sun seed
    sun <- which.max(membership_posibility)

    # Find circumsolar region
    ## get coordinates of every object
    rcells <- r
    rcells[] <- 1:ncell(r)
    .get_center <- function(x) {
      if (length(x) > 1) {
        xy <- terra::xyFromCell(r, x)
        x <- x[grDevices::chull(xy)]

        rad_z <- z[x] %>% .degree2radian()
        rad_a <- a[x] %>% .degree2radian()

        x3d <- sin(rad_z) * cos(rad_a)
        y3d <- sin(rad_z) * sin(rad_a)
        z3d <- cos(rad_z)

        x_mean <- mean(x3d[,])
        y_mean <- mean(y3d[,])
        z_mean <- mean(z3d[,])

        # to transform it to the unit vector in order to force it to be on the
        # sphere of radius one
        norm <- sqrt(x_mean^2 + y_mean^2 + z_mean^2)
        x_mean <- x_mean / norm
        y_mean <- y_mean / norm
        z_mean <- z_mean / norm

        # codify
        zenith <- acos(z_mean) %>% .radian2degree()
        azimuth <- atan2(y_mean, x_mean) %% (2 * pi) %>% .radian2degree()
        za <- c(zenith, azimuth)
      } else {
        za <- c(z[x][,], a[x][,])
      }
      round(za[1], 2) * 1e6 + round(za[2], 2)
    }
    za <-  compute_mean(rcells, labeled_m, .get_center, return = "vector")


    ## decodify
    zenith <- trunc(za/1e4)
    azimuth <- (za/1e4 - zenith) * 1e4
    zenith <- zenith/100
    za <- data.frame(zenith, azimuth) %>% .degree2radian()

    ## calc distance to sun seed
    d <- calc_spherical_distance(za[, 1], za[, 2], za[sun, 1], za[sun, 2])
    seg_labels <- compute_mean(labeled_m, labeled_m, return = "vector")

    ## classify circumsolar region based on distance
    i <- d > .degree2radian(circumsolar_angle)
    if (any(i)) {
      rcl <- data.frame(seg_labels[i], 0)
      m <- terra::classify(labeled_m, rcl)
    } else {
      m <- labeled_m
    }
    m <- m != 0

    # Calc coordinates of the circumsolar region
    za <-  compute_mean(rcells, m, .get_center, return = "vector")
    zenith <- trunc(za/1e4)
    azimuth <- (za/1e4 - zenith) * 1e4
    zenith <- zenith/100

    zenith_azimuth <- c(zenith, azimuth)
  } else {
    r[!bin] <- 0
    pan <- fisheye_to_pano(r, z, a,
                           function(x) quantile(x, 0.99, na.rm = TRUE), 1)
    zenith_azimuth <- terra::rowColFromCell(pan, which.max(pan[,])) %>%
      as.numeric()
  }
  sun_angles <- c(zenith_azimuth[1], zenith_azimuth[2])
  names(sun_angles) <- c("zen_deg", "azi_deg")
  sun_angles
}
