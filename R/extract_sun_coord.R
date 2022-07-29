#' Extract sun coordinates
#'
#' Extract the sun coordinates for CIE sky model fitting.
#'
#' This function uses an object-based image analyze theoretical framework. The
#' segmentation are given by \code{g} and \code{bin}. For every cell of
#' \code{g}, the pixels from \code{r} that are equal to one on \code{bin} are
#' selected, and its maximum value is calculated. Then, the 95th percentile of
#' these maximum values is computed and it is used to filter out cells below
#' that threshold; i.e, only the cells with at least one extremely bright sky
#' pixel is selected.
#'
#' The selected cells are grouped based on adjacency, and new bigger segments
#' are created from the groups. The degree of membership to the class \emph{Sun}
#' is calculated for every new segment by computing the number of cells that
#' constitute the segment and its mean digital number (values taken from
#' \code{r}). In other words, the brightest and largest segments are the ones
#' that score higher. The one with the highest score is selected as the
#' \emph{sun seed}.
#'
#' The angular distance from the sun seed to every other segments are computed,
#' and only the segments not farther than \code{max_angular_dist} are classified
#' as part of the sun corona. A multi-part segment is created by merging the
#' sun-corona segments and, finally, the center of its bounding box is
#' considered as the sun location.
#'
#' @inheritParams extract_sky_points
#' @inheritParams sky_grid_segmentation
#' @param max_angular_dist Numeric vector of length one. Angle in degree to
#'   establish the maximum size of the sun corona.
#'
#' @return Object of class \emph{list} with two numeric vectors of length two
#'   named \emph{row_col} and \emph{zenith_azimuth}. The former is raster
#'   coordinates of the solar disk (row and column), and the other is angular
#'   coordinates (zenith and azimuth angles in degrees).
#'
#' @family HSP Functions
#' @seealso fit_cie_sky_model
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' sky_blue_sample <- crop(caim, ext(610,643,760,806))
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' ecaim <- enhance_caim(caim, !is.na(z), sky_blue, gamma = 2.2)
#' bin <- apply_thr(ecaim, 0.75)
#' g <- sky_grid_segmentation(z, a, 10)
#' r <- gbc(caim$Blue*255)
#' sun_coord <- extract_sun_coord(r, z, a, bin, g, max_angular_dist = 30)
#' xy <- cellFromRowCol(z, sun_coord$row_col[1], sun_coord$row_col[2]) %>%
#'   xyFromCell(z, .)
#' plot(r)
#' plot(vect(xy), add = TRUE, col = 2)
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
           .calc_angular_distance(za[sun, 1], za[sun, 2], za[i, 1], za[i, 2]))
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
