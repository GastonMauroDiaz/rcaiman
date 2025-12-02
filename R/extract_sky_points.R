#' Extract sky points
#'
#' @description
#' Sample representative sky pixels for use in model fitting or interpolation.
#'
#' @details
#' Two sampling strategies are provided:
#' \describe{
#'   \item{`"per_cell"`}{select one sky point per cell of a segmentation (`seg`)
#'     as the brightest pixel marked `TRUE` in `bin`, provided the cell’s white
#'     pixel count exceeds one fourth of the mean across valid cells.}
#'   \item{`"local_max"`}{detect local maxima within a fixed \eqn{9 \times 9}
#'     window, restricted to pixels marked `TRUE` in `bin`, after removing
#'     patches of connected `TRUE` pixels that are implausible based on fixed
#'     area/size thresholds. Each detected maximum is taken as a sky point.
#'     }
#' }
#' Use `"per_cell"` to promote an even, representative spatial distribution (good
#' for model fitting), and `"local_max"` to be exhaustive for interpolation.
#'
#' @param r numeric [terra::SpatRaster-class] of one layer. Typically the blue
#'   band of a canopy image.
#' @param bin logical [terra::SpatRaster-class] of one layer. Binary image where
#'   `TRUE` marks candidate sky pixels. Typically the output of
#'   [binarize_with_thr()].
#' @param seg single-layer [terra::SpatRaster-class]. Segmentation map of r,
#'   typically created with functions such as [equalarea_segmentation()], [skygrid_segmentation()],
#'   [ring_segmentation()] or [sector_segmentation()], but any raster with
#'   integer segment labels is accepted. Ignored when `method = "local_max"`.
#' @param dist_to_black numeric vector of length one or `NULL`. Minimum distance
#'   (pixels) to the nearest black pixel for a candidate sky pixel to be valid.
#'   If `NULL`, no distance constraint is applied. Ignored with `method = local_max`.
#' @param method character vector of length one. Sampling method; either
#'   `"per_cell"` (default) or `"local_max"`.
#'
#' @return `data.frame` with columns `row` and `col`.
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
#' bin <- binarize_by_region(r, ring_segmentation(z, 15), "thr_isodata") &
#'   select_sky_region(z, 0, 88)
#'
#' seg <- skygrid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, seg,
#'                                  dist_to_black = 3)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' }
extract_sky_points <- function(r, bin, seg, dist_to_black = 3, method = "per_cell") {
  .this_requires_EBImage()

  .assert_single_layer(r)
  .assert_logical_mask(bin)
  .assert_same_geom(r, bin)
  .check_vector(dist_to_black, "integerish", 1, allow_null = TRUE, sign = "positive")
  .assert_choice(method, c("per_cell", "local_max"))

  if (method == "per_cell") {
    if (!is.null(dist_to_black)) {
      bin2 <- grow_black(bin, dist_to_black)
    } else {
      bin2 <- bin
    }
    .assert_single_layer(seg)
    .assert_same_geom(seg, r)

    # per-cell approach
    ## remove cells with not enough white pixels
    nwp <- extract_feature(bin, seg, sum, return = "vector")
    mean_nwp <- mean(nwp[nwp != 0])
    nwp <- extract_feature(bin, seg, sum, return = "raster")
    nwp <- nwp > mean_nwp/4
    seg[!nwp] <- 0

    no_col <- no_row <- bin
    terra::values(no_col) <- .col(dim(bin)[1:2])
    terra::values(no_row) <- .row(dim(bin)[1:2])

    ds <- data.frame(col = no_col[bin2],
                     row = no_row[bin2],
                     seg = seg[bin2],
                     dn = r[bin2])
    names(ds) <- c("col", "row", "seg", "dn")
    if (nrow(ds) != 0) { #to avoid crashing when there is no white pixels
        i <- tapply(1:nrow(ds), ds$seg,
                      function(x) {
                        x[which.max(ds$dn[x])]
                      })
        i <- i[names(i) != 0]
        ds <- ds[i,]
    }
    sky_points <- ds[, c("row", "col")]

  } else {
    # local-maximum approach

    bwlabels <- EBImage::bwlabel(as.array(bin))
    shape <- EBImage::computeFeatures.shape(bwlabels)
    bwlabels <- terra::setValues(bin, bwlabels)
    shape <- terra::subst(bwlabels, 1:nrow(shape), shape)

    ## Remove large gaps and artifacts
    shape[shape$s.area > 200 | shape$s.area == 0] <- NA

    ## Use effective circular area (ECA)
    eca <- shape$s.area / (shape$s.radius.sd + 1)
    bin[eca < 9] <- 0

    ## find local maximum
    lmax <- terra::focal(r*bin, 9, "max")
    lmax <- (lmax == r) * bin
    lmax[is.na(lmax)] <- 0
    i <- lmax[] %>% as.numeric() %>% as.logical()

    ## Create points
    sky_points <- rowColFromCell(lmax, cells(lmax)[i]) %>% as.data.frame()
    colnames(sky_points) <- c("row", "col")
  }

  sky_points
}
