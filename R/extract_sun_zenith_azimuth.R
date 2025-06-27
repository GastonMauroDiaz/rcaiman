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
#' @param g [SpatRaster-class] built with [sky_grid_segmentation()] or
#'   [chessboard()], or `NULL`. Thought this arguments the user can select
#'   between and object-based method when an actual sky grid is provided, and a
#'   pixel-based method, when `NULL` is provided.
#' @param chi_max_sun Numeric vector of length one. Specifies the maximum
#'   expected size of the circumsolar region in degrees.
#'
#' @return Numeric vector of length two, where the first element is the solar
#'   zenith angle and the second is the solar azimuth angle, both expressed in
#'   degrees.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' com <- compute_complementary_gradients(caim)
#' chroma <- max(com$blue_yellow, com$cyan_red)
#' bin <- apply_thr(chroma, thr_isodata(chroma[!is.na(chroma)]))
#' bin <- bin & apply_thr(com$blue_yellow, -0.2)
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_zenith_azimuth <- extract_sun_zenith_azimuth(r, z, a, bin, g,
#'                                                  chi_max_sun = 30)
#' row_col <- row_col_from_zenith_azimuth(z, a,
#'                                        sun_zenith_azimuth[1],
#'                                        sun_zenith_azimuth[2])
#' points(row_col[1,2], nrow(caim) - row_col[1,1], col = 3, pch = 10)
#' }
extract_sun_zenith_azimuth <- function(r, z, a, bin, g,
                                       chi_max_sun = 30) {
  .this_requires_EBImage()
  .check_if_r_z_and_a_are_ok(r, z, a)
  .is_single_layer_raster(bin, "bin")
  terra::compareGeom(bin, r)
  stopifnot(length(chi_max_sun) == 1)

  if (!is.null(g)) {
    .is_single_layer_raster(g, "g")
    terra::compareGeom(g, r)

    # Select cells with at least one extremely bright sky pixel
    g[!bin] <- NA
    r <- extract_feature(r, g, function(x) quantile(x, 0.99, na.rm = TRUE))
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
    size <- extract_feature(g, labeled_m, .fun, return_raster = FALSE) %>%
      normalize_minmax()
    dn <- extract_feature(r, labeled_m, median, return_raster = FALSE) %>%
      normalize_minmax()
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
    i <- d > .degree2radian(chi_max_sun)
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
    r[!bin] <- 0
    # r <- terra::focal(r, 3, mean)
    pan <- fisheye_to_pano(r, z, a,
                           function(x) quantile(x, 0.99, na.rm = TRUE), 1)
    zenith_azimuth <- terra::rowColFromCell(pan, which.max(pan[,])) %>%
      as.numeric()
  }
  zenith_azimuth
}
