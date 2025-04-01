#' Extract sky points
#'
#' Extract sky points for model fitting
#'
#' This function will automatically sample sky pixels from the sky regions
#' delimited by `bin`. The density and distribution of the sampling points is
#' controlled by the arguments `g`, `dist_to_black`, and `min_raster_dist`.
#'
#' As the first step, sky pixels from `r` are evaluated to find the pixel with
#' maximum digital value per cell of the `g` argument. The
#' `dist_to_black` argument allows users to establish a buffer zone for `bin`,
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
#' @param dist_to_black,min_raster_dist Numeric vector of length one or `NULL`.
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
#' sky_blue <- polarLAB(50, 17, 293)
#' ecaim <- enhance_caim(caim, m, sky_blue = sky_blue)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g,
#'                                  dist_to_black = 3,
#'                                  min_raster_dist = 10)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' }
extract_sky_points <- function(r, bin, g,
                               dist_to_black = 3,
                               min_raster_dist = 3) {

  .this_requires_EBImage()
  .is_single_layer_raster(r)
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  if (!is.null(g)) {
    .is_single_layer_raster(g, "g")
    terra::compareGeom(r, g)
  }
  terra::compareGeom(r, bin)
  if (!is.null(dist_to_black)) stopifnot(length(dist_to_black) == 1)
  if (!is.null(min_raster_dist)) stopifnot(length(min_raster_dist) == 1)

  if (!is.null(dist_to_black)) {
    stopifnot(.is_whole(dist_to_black))
    kern <- EBImage::makeBrush(dist_to_black, "box")
    bin2 <- bin
    bin2 <- EBImage::erode(as.array(bin2), kern) %>%
      terra::setValues(bin2, .)
    bin2[is.na(bin2)] <- 0
    bin2 <- as.logical(bin2)
  } else {
    bin2 <- bin
  }

  if (!is.null(g)) {
    # sky-grid approach

    ## remove cells with not enough white pixels
    nwp <- extract_feature(bin, g, sum, return_raster = FALSE)
    mean_nwp <- mean(nwp[nwp != 0])
    nwp <- extract_feature(bin, g, sum, return_raster = TRUE)
    nwp <- nwp > mean_nwp/4
    g[!nwp] <- 0

    no_col <- no_row <- bin
    terra::values(no_col) <- .col(dim(bin)[1:2])
    terra::values(no_row) <- .row(dim(bin)[1:2])

    ds <- data.frame(col = no_col[bin2],
                     row = no_row[bin2],
                     g = g[bin2],
                     dn = r[bin2])
    names(ds) <- c("col", "row", "g", "dn")
    if (nrow(ds) != 0) { #to avoid crashing when there is no white pixels
        i <- tapply(1:nrow(ds), ds$g,
                      function(x) {
                        x[which.max(ds$dn[x])]
                      })
        i <- i[names(i) != 0]
        ds <- ds[i,]
        if (!is.null(min_raster_dist) & nrow(ds) > 1) {
          ds <- .filter(ds, c("col", "row"), min_raster_dist)
        }
    }
    sky_points <- ds[,c("row", "col")]

  } else {
    # local-maximum approach

    bwlabels <- EBImage::bwlabel(as.array(bin2))
    shape <- EBImage::computeFeatures.shape(bwlabels)
    bwlabels <- terra::setValues(bin2, bwlabels)
    shape <- terra::subst(bwlabels, 1:nrow(shape), shape)

    ## Remove large gaps and artifacts
    shape[shape$s.area > 200 | shape$s.area == 0] <- NA

    ## Use effective circular area (ECA)
    eca <- shape$s.area / (shape$s.radius.sd + 1)
    bin2[eca < 9] <- 0

    ## find local maximum
    lmax <- terra::focal(r*bin2, 9, "max")
    lmax <- (lmax == r) * bin2
    lmax[is.na(lmax)] <- 0
    i <- lmax[] %>% as.numeric() %>% as.logical()

    ## Create points
    sky_points <- rowColFromCell(lmax, cells(lmax)[i]) %>% as.data.frame()
    colnames(sky_points) <- c("row", "col")
  }
  sky_points
}
