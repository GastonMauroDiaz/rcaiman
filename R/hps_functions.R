#' Extract sky marks
#'
#' Extract sky marks and create a special file to interface with HSP software.
#'
#' This function is part of a workflow that connects the \code{\link{mblt}}
#' algorithm with the HSP software package \insertCite{Lang2013}{rcaiman}. This
#' workflow aims to fully automatize the single camera linear ratio method
#' presented by \insertCite{Lang2010;textual}{rcaiman}.
#'
#' The \code{bin} argument should be a binarization from the \code{\link{mblt}}
#' method set with a high \code{w}. The latter is for obtaining pure sky pixels.
#' Setting \code{w} to maximum is not recommended because high frequency of
#' omission errors are expected with that setting.
#'
#' The pure sky pixels masked by \code{bin} will be automatically sampled by the
#' function. The density and distribution of the sampling points is controlled
#' by the arguments \code{g}, \code{dist_to_plant}, \code{angular_dist},
#' \code{min_raster_dist}, and \code{fuzzy_logic}.
#'
#' As first step, the maximum pure sky pixel is extracted from every cell of the
#' sky grid \code{g}. Then, if \code{dist_to_plant} is other than \code{NULL},
#' the distance to the non-pure sky pixels is computed for every pure sky
#' pixels, and them are filtrated based on \code{dist_to_plant}. This allows the
#' user to indicate the minimum gap size involved in the calculations. The
#' following steps allows to decrease the number of points while maintaining a
#' good distribution.
#'
#' If \code{fuzzy_logic} is set to \code{TRUE}, the parameter
#' \code{angular_dist} is used to compute how much the user considered too close
#' in the polar space. Then, starting with a given point, fuzzy logic is used to
#' filter out the darkest and closest points. If \code{fuzzy_logic} is set to
#' \code{FALSE}, the filter only consider distance and \code{angular_dist}
#' becomes a threshold.
#'
#' The \code{min_raster_dist} argument is a threshold and can be used to
#' indicate the minimum distance in the raster space that is required between
#' points.
#'
#' @param r \linkS4class{RasterLayer}. Greyscale image hemispherical image.
#' @param bin \linkS4class{RasterLayer}. The result of binarizing the \code{r}
#'   argument.
#' @param g \linkS4class{RasterLayer}. The result of a call to
#'   \code{\link{sky_grid_segmentation}} taking into account the camera, lens,
#'   and pre-processing involved in obtaining the \code{r} argument.
#' @param path_to_HPS_project Character vector of length one. Path to the HPS
#'   project folder.
#' @param img_name Character vector of length one. The name of the file
#'   associated to the \code{r} argument.
#' @param dist_to_plant Numeric vector of length one.
#' @param angular_dist Numeric vector of length one. Distance in degrees.
#' @param min_raster_dist Numeric vector of length one. Distance in pixels.
#' @param fuzzy_logic Logical vector of length one. Default is \code{TRUE}.
#'
#' @references
#' \insertRef{Lang2010}{rcaiman}
#'
#' \insertRef{Lang2013}{rcaiman}
#'
#' @seealso extract_sun_mark
#'
#' @return None. A file will be written in the HSP project folder.
#' @export
extract_sky_marks <- function(r, bin, g, path_to_HPS_project, img_name,
                              dist_to_plant = 3,
                              angular_dist = 5,
                              min_raster_dist = 9,
                              fuzzy_logic = TRUE ) {

  .is_integerish <- function(x) x == round(x)
  stopifnot(.is_integerish(dist_to_plant))


  # remove the pixels with NA neighbors because HSP extract with 3x3 kernel
  NA_count <- focal(!bin, matrix(1, 3, 3))

  no_col <- no_row <- bin
  no_col[] <- .col(dim(bin)[1:2])
  no_row[] <- .row(dim(bin)[1:2])


  # systematic sampling using a sky grid by taking the maximum from each cell
  if (!is.null(dist_to_plant)) {
    dist_to_plant_img <- NA_count == 0
    dist_to_plant_img[NA_count == 0] <- NA
    dist_to_plant_img <- distance(dist_to_plant_img)

    ds <- data.frame(col = no_col[dist_to_plant_img > dist_to_plant],
                     row = no_row[dist_to_plant_img > dist_to_plant],
                     g = g[dist_to_plant_img > dist_to_plant],
                     z = z[dist_to_plant_img > dist_to_plant],
                     a = a[dist_to_plant_img > dist_to_plant],
                     dn = r[dist_to_plant_img > dist_to_plant])

  } else {
    ds <- data.frame(col = no_col[bin],
                     row = no_row[bin],
                     g = g[bin],
                     z = z[bin],
                     a = a[bin],
                     dn = r[bin])
  }

  indices <- tapply(1:nrow(ds), ds$g, function(x) x[which.max(ds$dn[x])])
  ds <- ds[indices,]


  # filtering

  .filter_fuzzy <- function(ds, col_names, thr) {
    d <- as.matrix(dist(ds[, col_names]))
    indices <- c()
    i <- 0
    while (i < nrow(d)) {
      i <- i + 1
      indices <- c(indices, row.names(d)[i]) #include the point itself (p)
      x <- names(d[i, d[i,] <= thr])
      if (!is.null(x)) {
        dn_feature <- ds[x, "dn"]
        dn_feature <- normalize(dn_feature, min(dn_feature), max(dn_feature))
        dist_feature <- d[i,x]
        dist_feature <- normalize(dist_feature,
                                  min(dist_feature),
                                  max(dist_feature))
        dist_feature <- 1 - dist_feature^2
        feature <- dn_feature * dist_feature

        tmp <- x[feature > 0.75]
        if(length(tmp) != 0) x <- tmp
        # this exclude from future search the darkest points closer to p,
        # including itself
        rows2crop <- (1:nrow(d))[match(x, rownames(d))]
        cols2crop <- (1:ncol(d))[match(x, colnames(d))]
        d <- d[-rows2crop, -cols2crop]
      }
    }
    ds[indices,]
  }

  .filter <- function(ds, col_names, thr) {
    d <- as.matrix(dist(ds[, col_names]))
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

  if (!is.null(angular_dist)) {
    stopifnot(angular_dist >= 1)
    if (fuzzy_logic) {
      ds <- .filter_fuzzy(ds, c("z", "a"), angular_dist)
    } else {
      ds <- .filter(ds, c("z", "a"), angular_dist)
    }

  }
  if (!is.null(min_raster_dist)) {
    stopifnot(min_raster_dist >= 1)
    ds <- .filter(ds, c("col", "row"), min_raster_dist)
  }

  no <- nrow(ds)

  col.row_coordinates <- paste(ds$col, ds$row, "3", sep = ".")
  col.row_coordinates <- paste(col.row_coordinates, collapse = " ")

  ds <- c(no, col.row_coordinates)
  ds <- data.frame(ds)
  extension(img_name) <- ""
  write.table(ds, file.path(path_to_HPS_project,
                            "manipulate",
                            paste0(img_name, "_points.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}



#' Extract sun mark
#'
#' Extract sun mark and create a special file to interface with HSP software.
#'
#' @inheritParams extract_sky_marks
#'
#' @seealso extract_sky_marks
#' @export
extract_sun_mark <- function(r, bin, g, path_to_HPS_project, img_name) {
  r <- extract_feature(r, g, max)
  m <- r > quantile(r[], 0.99, na.rm = TRUE)
  m[is.na(g)] <- 0

  labeled_m <- EBImage::bwlabel(as.matrix(m))
  labeled_m <- raster(labeled_m)
  labeled_m[labeled_m == 0] <- NA
  browser()
  fun <- function(x) {
    x <- unique(x)
    length(x)
  }
  area <- extract_feature(g, labeled_m, fun, return_raster = FALSE) %>%
          normalize(., min(.), max(.))
  dn <- extract_feature(r, labeled_m, mean, return_raster = FALSE) %>%
        normalize(., min(.), max(.))
  membership_posibility <- area * dn
  sun <- which.max(membership_posibility)
  az
  features <- EBImage::computeFeatures.shape(labeled_m)
  features


  plot((labeled_m))


  no_col <- no_row <- r
  no_col[] <- .col(dim(r)[1:2])
  no_row[] <- .row(dim(r)[1:2])

  ds <- cbind(no_col[m], no_row[m])
  indexes <- chull(ds)
  p <- sp::SpatialPoints(ds[indexes,])
  p <- rgeos::gCentroid(p)
  sun <- paste(round(coordinates(p)), collapse = ".")
  extension(img_name) <- ""
  write.table(sun, file.path(path_to_HPS_project,
                             "manipulate",
                             paste0(img_name, "_sun.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")


}
