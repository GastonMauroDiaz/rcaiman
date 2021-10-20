#' Extract sky marks
#'
#' Extract sky marks for CIE sky model fitting
#'
#' The \code{bin} argument should be any binarized image that masked out pure
#' canopy (non-gap) pixels and mixed pixels, so to establish a region of
#' interest dominated by pure sky pixels (a.k.a., gap pixels).
#'
#' The \code{bin} argument can be obtained with the \code{\link{mblt}} method
#' set with a high \code{w}. The latter is for obtaining good candidates for the
#' pure sky pixel class. Setting \code{w} to maximum is not recommended because
#' that tend to produce high frequency of omission errors.
#'
#' This function will automatically sample in the sky region delimited by
#' \code{bin}. The density and distribution of the sampling points is controlled
#' by the arguments \code{g}, \code{dist_to_plant}, \code{raster_dist}, and
#' \code{fuzzy_logic}.
#'
#' As first step, the digital number under the class *Gap* --digital numbers
#' from \code{r} covered by pixels values equal to one on the \code{bin} layer--
#' are evaluated to extract its maximum value per cell of the sky grid \code{g}.
#' But, \code{dist_to_plant} allows users to establish a buffer zone before that
#' this extraction occur.
#'
#' The final step filter these local maximum values.If \code{fuzzy_logic} is set
#' to \code{TRUE}, the parameter \code{raster_dist} is used to compute how much
#' the user considered too near in the raster space. Then, starting with a given
#' point, fuzzy logic is used to filter out the nearest and darkest (lowest
#' digital number) points. If \code{fuzzy_logic} is set to \code{FALSE},
#' \code{angular_dist} becomes a distance threshold that ignore digital number.
#'
#'
#' @param r \linkS4class{RasterLayer}. Greyscale image hemispherical image.
#' @param bin \linkS4class{RasterLayer}. The result of binarizing the \code{r}
#'   argument.
#' @param g \linkS4class{RasterLayer}. The result of a call to
#'   \code{\link{sky_grid_segmentation}} taking into account the camera, lens,
#'   and pre-processing involved in obtaining the \code{r} argument.
#' @inheritParams sky_grid_segmentation
#' @param dist_to_plant Numeric vector of length one.
#' @param raster_dist Numeric vector of length one. Distance in pixels.
#' @param fuzzy_logic Logical vector of length one. Default is \code{TRUE}.
#'
#' @family hps functions
#'
#' @export
extract_sky_marks <- function(r, bin, g,
                              dist_to_plant = 3,
                              raster_dist = 9,
                              fuzzy_logic = TRUE ) {

  # remove the pixels with NA neighbors because HSP extract with 3x3 kernel
  NA_count <- focal(!bin, matrix(1, 3, 3))

  no_col <- no_row <- bin
  no_col[] <- .col(dim(bin)[1:2])
  no_row[] <- .row(dim(bin)[1:2])


  # systematic sampling using a sky grid by taking the maximum from each cell
  if (!is.null(dist_to_plant)) {
    stopifnot(.is_integerish(dist_to_plant))
    .this_requires_EBImage()
    kern <- EBImage::makeBrush(dist_to_plant, "box")
    dist_to_plant_img <- NA_count == 0
    dist_to_plant_img <- EBImage::erode(as.matrix(dist_to_plant_img), kern) %>%
                         setValues(dist_to_plant_img, .)

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

  if (!is.null(raster_dist)) {
    stopifnot(raster_dist >= 1)
    if (fuzzy_logic) {
      ds <- .filter_fuzzy(ds, c("col", "row"), raster_dist)
    } else {
      ds <- .filter(ds, c("col", "row"), raster_dist)
    }

  }

  ds[,-3]

}

#' Write sky marks
#'
#' Create a special file to interface with HSP software.
#'
#' This function is part of a workflow that connects the mblt algorithm with the
#' HSP software package \insertCite{Lang2013}{rcaiman}.
#'
#' This function was designed to be called after \code{\link{extract_sky_marks}}
#' in a workflow that connects the mblt algorithm with the HSP software package.
#' I such a workflow, the \code{\link{extract_sky_marks}} will use an image
#' pre-processed by the HSP software as its \code{r} argument. The
#' \code{img_name} argument of \code{write_sky_marks()} should be the name of
#' the file associated to the aforementioned \code{r} argument.
#'
#'
#' @param x Object from the class data.frame. The result of a calling to
#'   \code{\link{extract_sky_marks}}.
#' @param path_to_HSP_project Character vector of length one. Path to the HSP
#'   project folder.
#' @param img_name Character vector of length one. See details.
#'
#' @family hsp functions
#'
#' @references \insertRef{Lang2013}{rcaiman}
#'
#' @return None. A file will be written in the HSP project folder.
#' @export
#'
write_sky_marks <- function(x, path_to_HSP_project, img_name) {
  no <- nrow(ds)

  col.row_coordinates <- paste(ds$col, ds$row, "3", sep = ".")
  col.row_coordinates <- paste(col.row_coordinates, collapse = " ")

  ds <- c(no, col.row_coordinates)
  ds <- data.frame(ds)
  extension(img_name) <- ""
  write.table(ds, file.path(path_to_HSP_project,
                            "manipulate",
                            paste0(img_name, "_points.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}



#' Extract sun mark
#'
#' Extract the sun mark for CIE sky model fitting
#'
#' The \code{bin} argument should be the same than for
#' \code{\link{write_sky_marks}}
#'
#' @inheritParams extract_sky_marks
#'
#' @family hsp functions
#' @export
extract_sun_mark <- function(r, bin, g, path_to_HSP_project, img_name) {
  .this_requires_EBImage()

  r <- extract_feature(r, g, max)
  m <- r > quantile(r[], 0.99, na.rm = TRUE)
  m[is.na(g)] <- 0

  labeled_m <- EBImage::bwlabel(as.matrix(m))
  labeled_m <- raster(labeled_m)
  labeled_m[labeled_m == 0] <- NA

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

  angle_res <- 360 / round(.get_max(g) / 1000)
  a <- round(g / 1000) * angle_res
  azimuth <- extract_feature(a, labeled_m, mean, return_raster = FALSE)

  indices <- as.matrix(dist(azimuth))[sun,] > 90

  rcl <- data.frame(unique(labeled_m), unique(labeled_m))
  rcl[indices, 2] <- NA
  m <- reclassify(labeled_m, rcl)
  m <- !is.na(m)
  extent(m) <- extent(r)
  projection(m) <- NA

  no_col <- no_row <- r
  no_col[] <- .col(dim(r)[1:2])
  no_row[] <- .row(dim(r)[1:2])

  ds <- cbind(no_col[m], no_row[m])
  indexes <- chull(ds)
  p <- sp::SpatialPoints(ds[indexes,])
  p <- rgeos::gCentroid(p)
  col_row <- as.numeric(round(coordinates(p)))
  attr(col_row, "name") <- "col_row"
  col_row
}



#' Write sun mark
#'
#' Create a special file to interface with HSP software.
#'
#' Please, see the Details section of this function:
#' \code{\link{write_sky_marks}}.
#'
#' @param x Object from the class data.frame. The result of a calling to
#'   \code{\link{extract_sky_marks}}.
#' @inheritParams write_sky_marks
#'
#' @family hsp functions
#'
#' @return None. A file will be written in the HSP project folder.
#' @export
write_sun_mark <- function(x, path_to_HSP_project, img_name) {
  sun <- paste(x, collapse = ".")
  extension(img_name) <- ""
  write.table(sun, file.path(path_to_HSP_project,
                             "manipulate",
                             paste0(img_name, "_sun.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}
