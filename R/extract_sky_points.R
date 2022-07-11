#' Extract sky points
#'
#' Extract sky points for model fitting.
#'
#' The \code{bin} argument should be any binarized image that masked out pure
#' canopy (non-gap) pixels and most of the mixed pixels, so that it establish a
#' region of interest dominated by pure sky pixels (a.k.a., gap pixels). This
#' argument can be obtained with \code{\link{find_sky_pixels}} and
#' \code{\link{find_sky_pixels_nonnull_criteria}}, among other alternatives.
#'
#' This function will automatically sample sky pixels from the sky region
#' delimited by \code{bin}. The density and distribution of the sampling points
#' is controlled by the arguments \code{g}, \code{dist_to_plant}, and
#' \code{min_raster_dist}.
#'
#' As first step, the digital number under the class \emph{Gap} --digital
#' numbers from \code{r} covered by pixels values equal to one on the \code{bin}
#' layer-- are evaluated to extract its maximum value per cell of \code{g}. But,
#' \code{dist_to_plant} allows users to establish a buffer zone for \code{bin}.
#'
#' The final step filter these local maximum values as follow. The
#' \code{min_raster_dist} argument is a minimum distance threshold between
#' points that is applied in the raster space.
#'
#' Using code \code{NULL} as argument skip the filtering step in question.
#'
#' @inheritParams fit_trend_surface
#' @param g \linkS4class{SpatRaster}. The result of a call to
#'   \code{\link{sky_grid_segmentation}} taking into account the camera, lens,
#'   and pre-processing involved in obtaining the \code{r} argument.
#' @param dist_to_plant Numeric vector of length one or \code{NULL}.
#' @param min_raster_dist Numeric vector of length one or \code{NULL}. Distance
#'   in pixels.
#'
#' @family Tools functions
#'
#' @return An object of the class \emph{data.frame} with two columns named
#'   \emph{col}, \emph{row}.
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
#' sky_points <- extract_sky_points(r, bin, g)
#' cells <- cellFromRowCol(z, sky_points$row, sky_points$col)
#' hist(r[cells][,1])
#' xy <- xyFromCell(z, cells)
#' plot(r)
#' plot(vect(xy), add = TRUE, col = 2)
#' }
extract_sky_points <- function(r, bin, g,
                              dist_to_plant = 3,
                              min_raster_dist = 3) {
  .is_single_layer_raster(r)
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  .is_single_layer_raster(g, "g")
  terra::compareGeom(r, bin)
  terra::compareGeom(r, g)
  if (!is.null(dist_to_plant)) stopifnot(length(dist_to_plant) == 1)
  if (!is.null(min_raster_dist)) stopifnot(length(min_raster_dist) == 1)

  # remove the pixels with NA neighbors because HSP extract with 3x3 window
  NA_count <- terra::focal(!bin, w = 3, fun = "sum")

  no_col <- no_row <- bin
  terra::values(no_col) <- .col(dim(bin)[1:2])
  terra::values(no_row) <- .row(dim(bin)[1:2])


  # systematic sampling using a sky grid by taking the maximum from each cell
  if (!is.null(dist_to_plant)) {
    stopifnot(.is_integerish(dist_to_plant))
    .this_requires_EBImage()
    kern <- EBImage::makeBrush(dist_to_plant, "box")
    dist_to_plant_img <- NA_count == 0
    dist_to_plant_img <- EBImage::erode(as.matrix(dist_to_plant_img), kern) %>%
      terra::setValues(dist_to_plant_img, .)
    dist_to_plant_img[is.na(dist_to_plant_img)] <- 0
    dist_to_plant_img <- as.logical(dist_to_plant_img)

    ds <- data.frame(col = no_col[dist_to_plant_img],
                     row = no_row[dist_to_plant_img],
                     g = g[dist_to_plant_img],
                     dn = r[dist_to_plant_img])

  } else {
    ds <- data.frame(col = no_col[bin],
                     row = no_row[bin],
                     g = g[bin],
                     dn = r[bin])
  }

  indices <- tapply(1:nrow(ds), ds$g,
                    function(x) {
                      # browsser()
                      x[which.max(ds$dn[x])]

                      # indices <- ds$dn[x] >= quantile(ds$dn[x], 0.9)
                      # sample(x[indices], 1)

                      # dn <- ds$dn[x]
                      # indices <- dn >= quantile(dn, 0.9)
                      # i <- which.min(dn[indices])
                      # x[i]
                    })
  ds <- ds[indices,]


  # filtering
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
    }
    ds[indices,]
  }

  if (!is.null(min_raster_dist)) {
    ds <- .filter(ds, c("col", "row"), min_raster_dist)
  }
  ds[,c(2, 1)]
}
