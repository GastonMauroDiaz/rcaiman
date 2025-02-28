#' Extract sun coordinates
#'
#' Extract the sun coordinates for CIE sky model fitting.
#'
#' This function uses an object-based image analyze framework. The segmentation
#' is given by `g` and `bin`. For every cell of `g`, the pixels
#' from `r` that are equal to one on `bin` are selected and its
#' maximum value is calculated. Then, the 95th percentile of these maximum
#' values is computed and used to filter out cells below that threshold; i.e,
#' only the cells with at least one extremely bright sky pixel is selected.
#'
#' The selected cells are grouped based on adjacency, and new bigger segments
#' are created from these groups. The degree of membership to the class
#' *Sun* is calculated for every new segment by computing the number of
#' cells that constitute the segment and its mean digital number (values taken
#' from `r`). In other words, the  largest and brightest segments are the
#' ones that score higher. The one with the highest score is selected as the
#' *sun seed*.
#'
#' The angular distance from the sun seed to every other segments are computed,
#' and only the segments not farther than `max_angular_dist` are classified
#' as part of the sun corona. A multi-part segment is created by merging the
#' sun-corona segments and, finally, the center of its bounding box is
#' considered as the sun location.
#'
#' @inheritParams extract_sky_points
#' @inheritParams sky_grid_segmentation
#' @param max_angular_dist Numeric vector of length one. Angle in degrees to
#'   control the potential maximum size of the solar corona.
#'
#' @return Object of class *list* with two numeric vectors of length two
#'   named *row_col* and *zenith_azimuth*. The former is the raster
#'   coordinates of the solar disk (row and column), and the other is the
#'   angular coordinates (zenith and azimuth angles in degrees).
#'
#' @family Tool Functions
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
#' mx <- optim_normalize(caim, bin)
#' caim <- normalize(caim, 0, mx, TRUE)
#' plotRGB(caim*255)
#' sky_blue <- HSV(239, 0.85, 0.5)
#' ecaim <- enhance_caim(caim, m, sky_blue = sky_blue)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_coord <- extract_sun_coord(r, z, a, bin, g, max_angular_dist = 30)
#' points(sun_coord$row_col[2], nrow(caim) - sun_coord$row_col[1],
#'         col = 3, pch = 10)
#' }
extract_sun_coord <- function(r, z, a, bin, g,
                             max_angular_dist = 30) {
  .this_requires_EBImage()
  .check_if_r_z_and_a_are_ok(r, z, a)
  .is_single_layer_raster(bin, "bin")
  .is_single_layer_raster(g, "g")
  terra::compareGeom(bin, r)
  terra::compareGeom(g, r)
  stopifnot(length(max_angular_dist) == 1)

  g[!bin] <- NA
  r <- extract_feature(r, g, max)
  m <- r >= quantile(terra::values(r), 0.95, na.rm = TRUE)
  m[is.na(m)] <- 0

  labeled_m <- EBImage::bwlabel(as.array(m)[,,1])
  labeled_m <- terra::rast(labeled_m)
  terra::ext(labeled_m) <- terra::ext(r)
  terra::crs(labeled_m) <- terra::crs(r)

  fun <- function(x) {
    x <- unique(x)
    length(x)
  }
  area <- extract_feature(g, labeled_m, fun, return_raster = FALSE) %>%
    normalize()
  dn <- extract_feature(r, labeled_m, mean, return_raster = FALSE) %>%
    normalize()
  if (any(is.nan(dn))) dn[] <- 1
  if (any(is.nan(area))) area[] <- 1
  membership_posibility <- area * dn
  sun <- which.max(membership_posibility)

  fun <- function(x) mean(range(x))
  azimuth <- extract_feature(a, labeled_m, fun, return_raster = FALSE) %>%
    .degree2radian()
  zenith <- extract_feature(z, labeled_m, fun, return_raster = FALSE) %>%
    .degree2radian()
  za <- data.frame(zenith, azimuth)
  d <- c()
  for (i in 1:nrow(za)) {
    d <- c(d,
           .calc_spherical_distance(za[sun, 1], za[sun, 2], za[i, 1], za[i, 2]))
  }

  indices <- d > .degree2radian(max_angular_dist)
  rcl <- data.frame(names(zenith) %>% as.numeric(),
                    names(zenith) %>% as.numeric())
  rcl[indices, 2] <- 0
  m <- terra::classify(labeled_m, rcl)
  m <- m != 0

  no_col <- no_row <- r
  terra::values(no_col) <- .col(dim(r)[1:2])
  terra::values(no_row) <- .row(dim(r)[1:2])

  row_col <- data.frame(extract_feature(no_row, m, fun, return_raster = FALSE),
                        extract_feature(no_col, m, fun, return_raster = FALSE))
  row_col <- round(row_col) %>% as.numeric()
  sun_coord <- c(z[row_col[1], row_col[2]] %>% as.numeric(),
                 a[row_col[1], row_col[2]] %>% as.numeric())
  list(row_col = row_col, zenith_azimuth = round(sun_coord))
}
