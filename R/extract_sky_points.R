#' Extract sky points
#'
#' Extract sky points for model fitting
#'
#' This function will automatically sample sky pixels from the sky regions
#' delimited by `bin`. The density and distribution of the sampling points is
#' controlled by the arguments `g`, `dist_to_plant`, and `min_raster_dist`.
#'
#' As the first step, sky pixels from `r` are evaluated to find the pixel with
#' maximum digital value (local maximum) per cell of the `g` argument. The
#' `dist_to_plant` argument allows users to establish a buffer zone for `bin`,
#' meaning a size reduction of the original sky regions.
#'
#' The final step is filtering these local maximum values by evaluating the
#' Euclidean distances between them on the raster space. Any new point with a
#' distance from existing points minor than `min_raster_dist` is discarded. Cell
#' labels determine the order in which the points are evaluated.
#'
#' To skip a given filtering step, use code `NULL` as argument input. For
#' instance, `min_raster_dist = NULL` will return points omitting
#' the final step.
#'
#' @inheritParams fit_trend_surface
#' @param g [SpatRaster-class] built with [sky_grid_segmentation()] or
#'   [chessboard()].
#' @param dist_to_plant,min_raster_dist Numeric vector of length one or `NULL`.
#'
#'
#' @family Tool Functions
#' @seealso [fit_cie_sky_model()]
#'
#' @return An object of the class *data.frame* with two columns named
#'   *row* and *col*.
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
#' sky_points <- extract_sky_points(r, bin, g,
#'                                  dist_to_plant = 3,
#'                                  min_raster_dist = 10)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' }
extract_sky_points <- function(r, bin, g,
                              dist_to_plant = 3,
                              min_raster_dist = 3) {


  .filter <- function(ds, col_names, thr) {
    d <- as.matrix(stats::dist(ds[, col_names]))
    indices <- c()
    i <- 0
    while (i < nrow(d)) {
      i <- i + 1
      indices <- c(indices, row.names(d)[i]) #include the point itself (p)
      x <- names(d[i, d[i,] <= thr])
      if (!is.null(x)) {
        # this exclude from future search all the points near p,
        # including itself
        rows2crop <- (1:nrow(d))[match(x, rownames(d))]
        cols2crop <- (1:ncol(d))[match(x, colnames(d))]
        d <- d[-rows2crop, -cols2crop]
      }
      if (is.vector(d)) d <- matrix(d)
    }
    ds[indices,]
  }


  .is_single_layer_raster(r)
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  .is_single_layer_raster(g, "g")
  terra::compareGeom(r, bin)
  terra::compareGeom(r, g)
  if (!is.null(dist_to_plant)) stopifnot(length(dist_to_plant) == 1)
  if (!is.null(min_raster_dist)) stopifnot(length(min_raster_dist) == 1)

  #new feature
  # cv <- terra::focal(r, 3, sd) / terra::focal(r, 3, mean)
  # r <- normalize(-cv) * r

  #new feature
  nw <- extract_feature(bin, g, sum, return_raster = FALSE)
  tol <- sum(nw != 0) / 3
  size <- extract_feature(g, g, length, return_raster = FALSE) %>% min()
  w <- 1
  while (sum(nw > size) < tol & w > 0) {
    w <- w - 0.01
    size <- size * w
  }
  nw <- extract_feature(bin, g, sum, return_raster = TRUE)
  nw <- nw > size
  g[!nw] <- 0

  # remove the pixels with NA neighbors because HSP extract with 3x3 window
  # NA_count <- terra::focal(!bin, w = 3, fun = "sum")

  no_col <- no_row <- bin
  terra::values(no_col) <- .col(dim(bin)[1:2])
  terra::values(no_row) <- .row(dim(bin)[1:2])

  # systematic sampling using a sky grid by taking the maximum from each cell
  if (!is.null(dist_to_plant)) {
    stopifnot(.is_integerish(dist_to_plant))
    .this_requires_EBImage()
    kern <- EBImage::makeBrush(dist_to_plant, "box")
    # dist_to_plant_img <- NA_count == 0
    dist_to_plant_img <- bin
    dist_to_plant_img <- EBImage::erode(as.array(dist_to_plant_img), kern) %>%
      terra::setValues(dist_to_plant_img, .)
    dist_to_plant_img[is.na(dist_to_plant_img)] <- 0
    dist_to_plant_img <- as.logical(dist_to_plant_img)

    ds <- data.frame(col = no_col[dist_to_plant_img],
                     row = no_row[dist_to_plant_img],
                     g = g[dist_to_plant_img],
                     dn = r[dist_to_plant_img])
    names(ds) <- c("col", "row", "g", "dn")
  } else {
    ds <- data.frame(col = no_col[bin],
                     row = no_row[bin],
                     g = g[bin],
                     dn = r[bin])
    names(ds) <- c("col", "row", "g", "dn")
  }

  if (nrow(ds) != 0) {
      i <- tapply(1:nrow(ds), ds$g,
                    function(x) {
                      x[which.max(ds$dn[x])]
                    })
      ds <- ds[i,]
      if (!is.null(min_raster_dist) & nrow(ds) > 1) {
        ds <- .filter(ds, c("col", "row"), min_raster_dist)
      }
  }
  ds[,c(2, 1)]
}
