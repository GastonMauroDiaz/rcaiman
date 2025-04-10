#' Extract sun coordinates
#'
#' Extract the sun coordinates for CIE sky model fitting.
#'
#' This function uses an object-based image analyze framework. The segmentation
#' is given by `g` and `bin`. For every cell of `g`, the pixels from `r` that
#' are equal to one on `bin` are selected and its maximum value is calculated.
#' Then, the 95th percentile of these maximum values is computed and used to
#' filter out cells below that threshold; i.e, only the cells with at least one
#' extremely bright sky pixel is selected.
#'
#' The selected cells are grouped based on adjacency, and new bigger segments
#' are created from these groups. The degree of membership to the class
#' *Sun* is calculated for every new segment by computing the number of
#' cells that constitute the segment and its mean digital number (values taken
#' from `r`). In other words, the  largest and brightest segments are the ones
#' that score higher. The one with the highest score is selected as the
#' *sun seed*.
#'
#' The angular distance from the sun seed to every other segments are computed,
#' and only the segments not farther than `max_angular_dist` are classified as
#' part of the circumsolar region. A multi-part segment is created by merging
#' the circumsolar segments and, finally, the center of its bounding box is
#' considered as the sun location.
#'
#' @inheritParams extract_sky_points
#' @inheritParams sky_grid_segmentation
#' @param max_angular_dist Numeric vector of length one. Specifies the maximum
#'   expected size of the circumsolar region in degrees.
#'
#' @return Object of class *list* with two numeric vectors of length two named
#'   *row_col* and *zenith_azimuth*. The former is the raster coordinates of the
#'   solar disk (row and column), and the other is the angular coordinates
#'   (zenith and azimuth angles in degrees).
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
#' caim <- normalize_minmax(caim, 0, mx, TRUE)
#' plotRGB(caim*255)
#' sky_blue <- polarLAB(50, 17, 293)
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
  dn <- extract_feature(r, labeled_m, mean, return_raster = FALSE) %>%
    normalize_minmax()
  if (any(is.nan(dn))) dn[] <- 1
  if (any(is.nan(size))) size[] <- 1
  membership_posibility <- size * dn

  # Find sun seed
  sun <- which.max(membership_posibility)

  # Find circumsolar region
  ## get coordinates of every object
  rcells <- r
  rcells[] <- 1:ncell(r)
  .get_center <- function(x) {
    xy <- terra::xyFromCell(r, x)
    if (nrow(xy) > 1) {
      xy <- xy[grDevices::chull(xy),]
      v <- terra::vect(list(xy), type = "polygon", crs = crs(r))
      v <- terra::centroids(v)
      xy <- terra::crds(v)
    }
    terra::cellFromXY(r, xy)
  }
  cells <-  extract_feature(rcells, labeled_m, .get_center, return_raster = FALSE)
  za <- data.frame(zenith = z[cells], azimuth = a[cells]) %>% .degree2radian()

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
  cells <-  extract_feature(rcells, m, .get_center, return_raster = FALSE)
  row_col <- terra::rowColFromCell(r, cells)

  row_col <- round(row_col) %>% as.numeric()
  sun_coord <- c(z[row_col[1], row_col[2]] %>% as.numeric(),
                 a[row_col[1], row_col[2]] %>% as.numeric())


  list(row_col = row_col,
       zenith_azimuth = sun_coord)
}
