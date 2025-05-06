#' Extract sun angular coordinates
#'
#' Extract the sun coordinates for CIE sky model fitting.
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly presenting and
#' testing this pipeline is under preparation.
#'
#' @inheritParams extract_sky_points
#' @inheritParams sky_grid_segmentation
#' @param max_angular_dist Numeric vector of length one. Specifies the maximum
#'   expected size of the circumsolar region in degrees.
#' @param method description
#'
#' @return Numeric vector of length two, where the first element is the
#'   solar zenith angle and the second is the solar azimuth angle, both
#'   expressed in degrees.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' mx <- optim_max(caim, bin)
#' caim <- normalize_minmax(caim, 0, mx, TRUE)
#' plotRGB(caim*255)
#' sky_blue <- polarLAB(50, 17, 293)
#' ecaim <- enhance_caim(caim, m, sky_blue = sky_blue)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_zenith_azimuth <- extract_sun_zenith_azimuth(r, z, a, bin, g,
#'                                                  max_angular_dist = 30)
#' row_col <- row_col_from_zenith_azimuth(z, a,
#'                                        sun_zenith_azimuth[1],
#'                                        sun_zenith_azimuth[2])
#' points(row_col[1,2], nrow(caim) - row_col[1,1], col = 3, pch = 10)
#' }
extract_sun_zenith_azimuth <- function(r, z, a, bin, g,
                                       max_angular_dist = 30,
                                       method = "object-based") {
  .this_requires_EBImage()
  .check_if_r_z_and_a_are_ok(r, z, a)
  .is_single_layer_raster(bin, "bin")
  .is_single_layer_raster(g, "g")
  terra::compareGeom(bin, r)
  terra::compareGeom(g, r)
  stopifnot(length(max_angular_dist) == 1)

  stopifnot(method %in% c("object-based", "pixel-based"))

  if (method == "object-based") {
    # Select cells with at least one extremely bright sky pixel
    g[!bin] <- NA
    r <- extract_feature(r, g, max)
    m <- r >= quantile(terra::values(r), 0.95, na.rm = TRUE)
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
    size <- extract_feature(g, labeled_m, .fun, return_raster = FALSE) %>%
      normalize_minmax()
    dn <- extract_feature(r, labeled_m, median, return_raster = FALSE) %>%
      normalize_minmax()
    if (any(is.nan(dn))) dn[] <- 1
    if (any(is.nan(size))) size[] <- 1
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

        zenith <- acos(z_mean) %>% .radian2degree()
        azimuth <- atan2(y_mean, x_mean) %% (2 * pi) %>% .radian2degree()
        za <- c(zenith, azimuth)
      } else {
        za <- c(z[x][,], a[x][,])
      }
      round(za[1], 2) * 1e6 + round(za[2], 2)
    }
    za <-  extract_feature(rcells, labeled_m, .get_center, return_raster = FALSE)
    zenith <- trunc(za/1e4)
    azimuth <- (za/1e4 - zenith) * 1e4
    zenith <- zenith/100
    za <- data.frame(zenith, azimuth) %>% .degree2radian()


    ## calc distance to sun seed
    d <- calc_spherical_distance(za[, 1], za[, 2], za[sun, 1], za[sun, 2])
    seg_labels <- extract_feature(labeled_m, labeled_m, return_raster = FALSE)

    ## classify circumsolar region based on distance
    i <- d > .degree2radian(max_angular_dist)
    if (any(i)) {
      rcl <- data.frame(seg_labels[i], 0)
      m <- terra::classify(labeled_m, rcl)
    } else {
      m <- labeled_m
    }
    m <- m != 0

    # Calc coordinates of the circumsolar region
    za <-  extract_feature(rcells, m, .get_center, return_raster = FALSE)
    zenith <- trunc(za/1e4)
    azimuth <- (za/1e4 - zenith) * 1e4
    zenith <- zenith/100

    zenith_azimuth <- c(zenith, azimuth)
  } else {
    mat <- terra::focalMat(r, terra::ncol(r)/50, "circle")
    r_smoothed <- terra::focal(r, mat)
    m_sun <- (r_smoothed - max(r_smoothed[], na.rm = TRUE)) == 0

    zenith_azimuth <- c(z[m_sun][1,], a[m_sun][1,]) %>% unname()
  }
  zenith_azimuth
}
