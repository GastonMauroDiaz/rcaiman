#' Extract sky marks
#'
#' Extract sky marks for CIE sky model fitting
#'
#' The \code{bin} argument should be any binarized image that masked out pure
#' canopy (non-gap) pixels and most of the mixed pixels, so to establish a
#' region of interest dominated by pure sky pixels (a.k.a., gap pixels). The
#' \code{bin} argument can be obtained with the \code{\link{ootb_mblt}}.
#'
#' This function will automatically sample in the sky region delimited by
#' \code{bin}. The density and distribution of the sampling points is controlled
#' by the arguments \code{g}, \code{dist_to_plant}, and \code{min_raster_dist}.
#'
#' As first step, the digital number under the class *Gap* --digital numbers
#' from \code{r} covered by pixels values equal to one on the \code{bin} layer--
#' are evaluated to extract its maximum value per cell of the sky grid \code{g}.
#' But, \code{dist_to_plant} allows users to establish a buffer zone which
#' modify \code{bin}.
#'
#' The final step filter these local maximum values.\code{min_raster_dist} is a
#' minimum distance threshold between points that is applied in the raster
#' space.
#'
#' @param r \linkS4class{RasterLayer}. Greyscale image hemispherical image.
#' @param bin \linkS4class{RasterLayer}. The result of binarizing the \code{r}
#'   argument.
#' @param g \linkS4class{RasterLayer}. The result of a call to
#'   \code{\link{sky_grid_segmentation}} taking into account the camera, lens,
#'   and pre-processing involved in obtaining the \code{r} argument.
#' @param dist_to_plant Numeric vector of length one or \code{NULL}.
#' @param min_raster_dist Numeric vector of length one or \code{NULL}. Distance
#'   in pixels.
#'
#' @family hps functions
#'
#' @export
#' @examples
#' \dontrun{
#' require(magrittr)
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' blue <- gbc(r$Blue)
#' bin <- ootb_mblt(blue, z, a)$bin
#' sky_marks <- extract_sky_marks(blue, bin, g,
#'                                min_raster_dist = 10)
#' xy <- cellFromRowCol(z, sky_marks$row, sky_marks$col) %>%  xyFromCell(z, .)
#' plot(blue)
#' plot(SpatialPoints(xy), add = TRUE, col = 2)
#' }
extract_sky_marks <- function(r, bin, g,
                              dist_to_plant = 3,
                              min_raster_dist = 3) {

  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(bin) == "RasterLayer")
  stopifnot(class(g) == "RasterLayer")

  if (!is.null(dist_to_plant)) stopifnot(length(dist_to_plant) == 1)
  if (!is.null(min_raster_dist)) stopifnot(length(min_raster_dist) == 1)

  # remove the pixels with NA neighbors because HSP extract with 3x3 window
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

#' Write sky marks
#'
#' Create a special file to interface with HSP software.
#'
#' This function is part of a workflow that connects the MBLT algorithm with the
#' HSP software package \insertCite{Lang2013}{rcaiman}.
#'
#' This function was designed to be called after \code{\link{extract_sky_marks}}
#' in a workflow that connects the MBLT algorithm with the HSP software package.
#' I such a workflow, the \code{\link{extract_sky_marks}} will use an image
#' pre-processed by the HSP software as its \code{r} argument. Those images are
#' stored as PGM files by HSP in the subfolder "manipulate" of the project
#' folder (which will be in turn a subfolder of "project*s*" folder). Those PGM
#' files can be read with \code{\link[raster]{raster}}. For instance: \code{r <-
#' raster("C:/Users/johndoe/Documents/HSP/Projects/my_project/manipulate/img_001.pgm")}.
#' Then, they will be ready to use as input after running \code{normalize(r,
#' min(r[]),  max(r[]))}.
#'
#' The \code{img_name} argument of \code{write_sky_marks()} should be the name
#' of the file associated to the aforementioned \code{r} argument.
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
  utils::write.table(ds, file.path(path_to_HSP_project,
                            "manipulate",
                            paste0(img_name, "_points.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}



#' Extract sun mark
#'
#' Extract the sun mark for CIE sky model fitting.
#'
#' This function use an object-based image analyze theoretical framework. The
#' segmentation are given by \code{g}, and \code{b} should be understood as an
#' extra layer. For every cell of \code{g}, the maximum is calculated from the
#' pixel values on \code{r} that have a value equal to one on \code{bin}. Then,
#' the quantile 0.95 is extracted from these maximum values and is used to
#' filter out cells below that threshold, i.e, only the ones with at least one
#' extremely bright sky pixel is keep.
#'
#' Selected cells are grouped into segments based on adjacency. The degree of
#' membership to the class _Sun_ is calculated for every segment by using linear
#' membership functions for the features brightness --digital number from
#' \code{r}-- and size -- number of cells that constitute the segment. In other
#' words, the brighteners and lagers segments are the ones that score higher.
#' The one with the highest score is selected as a sun seed.
#'
#' The angular distance from the sun seed to every other segments are computed,
#' and only the segments not farther than \code{max_angular_dist} are classified
#' as part of the sun corona. A multi-part segment is created by merging the
#' sun-corona segments, a convex hull is calculated, and the centroide of its
#' shape is returned as the sun location in raster coordinates (column and row).
#'
#' The \code{bin} argument should be the same than for
#' \code{\link{extract_sky_marks}}.
#'
#' @inheritParams extract_sky_marks
#' @inheritParams sky_grid_segmentation
#' @param max_angular_dist Numeric vector of length one. Angle in degree to
#'   establish the maximum size of the sun corona. See details. The former is
#'   raster coordinates of the solar disk (row and column), and the other is
#'   angular coordinates (zenith and azimuth angles in degrees).
#'
#' @return Object of class list with two elements names "row_col" and
#'   "zenith_azimuth", both are numeric vectors of lenght two.
#'
#' @family hsp functions
#' @export
extract_sun_mark <- function(r, bin, z, a, g,
                             max_angular_dist = 45) {
  .this_requires_EBImage()
  .check_if_r_z_and_a_are_ok(r, z, a)

  stopifnot(class(bin) == "RasterLayer")
  stopifnot(class(g) == "RasterLayer")
  stopifnot(length(max_angular_dist) == 1)

  r <- extract_feature(r, g, max)
  m <- r >= quantile(r[], 0.95, na.rm = TRUE)
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
  if (any(is.nan(dn))) dn[] <- 1
  if (any(is.nan(area))) area[] <- 1
  membership_posibility <- area * dn
  sun <- which.max(membership_posibility)

  azimuth <- extract_feature(a, labeled_m, mean, return_raster = FALSE) %>%
    degree2radian()
  zenith <- extract_feature(z, labeled_m, mean, return_raster = FALSE) %>%
    degree2radian()
  za <- data.frame(zenith, azimuth)

  ids_to_vs <- expand.grid(sun, 1:nrow(za))

  d <- c()
  for (i in 1:nrow(za)) {
    d <- c(d, .calc_angular_distance(za[sun, 1], za[sun, 2], za[i, 1], za[i, 2]))
  }

  indices <- d > degree2radian(max_angular_dist)

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
  xy <- cellFromRowCol(m, no_row[m], no_col[m]) %>% xyFromCell(m, .)
  indexes <- grDevices::chull(xy)
  p <- sp::SpatialPoints(xy[indexes, ])

  p <- rgeos::gCentroid(p)
  xy <- as.numeric(round(coordinates(p)))
  row_col <- c(rowFromY(m, xy[2]), colFromX(m, xy[1]))

  sun_coord <- (c(z[row_col[1], row_col[2]], a[row_col[1], row_col[2]]))
  list(row_col = row_col, zenith_azimuth = round(sun_coord))
}



#' Write sun mark
#'
#' Create a special file to interface with HSP software.
#'
#' Refer to the Details section of this function:
#' \code{\link{write_sky_marks}}.
#'
#' @param x Numeric vector of lenght two. Raster coordinates of the solar disk
#'   that can be obtained by calling to \code{\link{extract_sun_mark}}. *TIP*:
#'   if the output of \code{extrac_sun_mark()} is \code{x}, then you should
#'   provide to \code{write_sun_mark()} this: \code{x$row_col}. See also
#'   \code{\link{row_col_from_zenith_azimuth}}.
#' @inheritParams write_sky_marks
#'
#'
#' @family hsp functions
#'
#' @return None. A file will be written in the HSP project folder.
#' @export
write_sun_mark <- function(x, path_to_HSP_project, img_name) {
  sun <- paste(x[c(2,1)], collapse = ".")
  extension(img_name) <- ""
  utils::write.table(sun, file.path(path_to_HSP_project,
                                    "manipulate",
                                    paste0(img_name, "_sun.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}


#' Row and col numbers from zenith and azimuth angles
#'
#' @inheritParams ootb_mblt
#' @param za Numeric vector of length two. Zenith and azimuth angles in degrees.
#'
#' @export
#'
#' @family hps functions
#'
#' @examples
#' z <- zenith_image(1000, lens_coef = lens())
#' row_col <- row_col_from_zenith_azimuth(z, c(45, 270), lens())
row_col_from_zenith_azimuth <- function(r, za, lens_coef) {
  stopifnot(class(r) == "RasterLayer")
  stopifnot(ncol(r) == nrow(r))
  stopifnot(is.numeric(lens_coef))
  stopifnot(is.numeric(za))
  stopifnot(length(za) == 2)
  az <- rev(za)
  rr <- calc_relative_radius(az[2], lens_coef)
  pol <- data.frame(theta = az[1] * pi/180 + pi/2,
                    r = rr * 90 * pi/180,
                    z = 0)
  cart <- pracma::pol2cart(as.matrix(pol))
  p <- sp::SpatialPoints(matrix(cart[1:2], ncol = 2))
  e <- extent(r)
  extent(r) <- extent(-pi/2,pi/2,-pi/2,pi/2)
  ir <- rasterize(p, r)
  cellsFromExtent(r, extent(p)) %>% rowColFromCell(r, .) %>% as.numeric()
}


#' Zenith and azimuth angles from row and col numbers
#'
#' @inheritParams ootb_mblt
#' @param row_col Numeric vector of length two. Row and col numbers.
#'
#' @export
#'
#' @family hps functions
#' @examples
#' z <- zenith_image(1000, lens_coef = lens())
#' zenith_azimuth_from_row_col(z, c(501, 751), lens())
zenith_azimuth_from_row_col <- function(r, row_col, lens_coef) {
  stopifnot(class(r) == "RasterLayer")
  stopifnot(ncol(r) == nrow(r))
  stopifnot(is.numeric(lens_coef))
  stopifnot(is.numeric(row_col))
  stopifnot(length(row_col) == 2)

  #get azimuth
  e <- extent(r)
  extent(r) <- extent(-pi/2,pi/2,-pi/2,pi/2)
  xy <- cellFromRowCol(r, row_col[1], row_col[2]) %>% xyFromCell(r, .)
  tr <- pracma::cart2pol(as.numeric(xy))
  azimuth <- tr[1] - pi/2 * 180/pi
  if (azimuth < 0) azimuth <- 360 + azimuth
  #get relative radius
  rr <- tr[2] * 180/pi / 90
  #invert
  zs <- seq(0,150, 0.1)
  rrs <- calc_relative_radius(zs, lens_coef)
  z_from_rr <- suppressWarnings(splinefun(rrs, zs))
  zenith <- z_from_rr(rr)
  c(zenith, azimuth)
}


#' Read manual input
#'
#' Read manual input stored in an HSP project.
#'
#' Refer to the Details section of this function:
#' \code{\link{write_sky_marks}}.
#'
#' @inheritParams write_sky_marks
#' @family hps functions
#' @export
read_manual_input <- function(path_to_HSP_project, img_name) {
  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "settings", full.names = TRUE)
  file <- files[grep(img_name, files)]
  settings <- scan(file, "character")
  settings <- settings[c(9, 11:13)]
  settings <- data.frame(
    name = Map(function(x) x[1], strsplit(settings, "=")) %>% unlist(),
    value = Map(function(x) x[2], strsplit(settings, "=")) %>% unlist()
  )

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "sun", full.names = TRUE)
  file <- files[grep(img_name, files)]
  sun <- scan(file, "character")
  sun <- strsplit(sun, "\\.") %>% unlist() %>% as.numeric()
  sun_mark <- list()
  sun_mark$row_col <- rev(sun)

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "points", full.names = TRUE)
  file <- files[grep(img_name, files)]
  sky_marks <- scan(file, "character", skip = 1)
  sky_marks <- strsplit(sky_marks, "\\.") %>%
    unlist() %>%
    as.numeric() %>%
    matrix(., ncol = 3, byrow = TRUE) %>%
    as.data.frame(.)
  names(sky_marks) <- c("col", "row", "type" )

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "statistics", full.names = TRUE)
  file <- files[grep(img_name, files)]
  content <- scan(file, "character", skip = 1, sep = "\n")
  zenith_dn <- content[grep( "Zenith", content)]
  zenith_dn <- strsplit(zenith_dn, "=")[[1]][2] %>%
    sub(",", ".", .) %>% as.numeric()

  list(weight = settings[1,2] %>% as.numeric(),
       max_points = settings[2,2] %>% as.numeric(),
       angle = settings[3,2] %>% as.numeric(),
       point_radius = settings[4,2] %>% as.numeric(),
       sun_mark = sun_mark,
       sky_marks = sky_marks,
       zenith_dn = zenith_dn)
}

#' Read optimized sky coefficients
#'
#' Read optimized CIE sky coefficients stored in an HSP project.
#'
#' Refer to the Details section of this function:
#' \code{\link{write_sky_marks}}.
#'
#' @inheritParams write_sky_marks
#' @family hps functions
#' @seealso cie_sky_model_raster
#' @export
read_opt_sky_coef <- function(path_to_HSP_project, img_name) {
  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "opt-parameters", full.names = TRUE)
  file <- files[grep(img_name, files)]
  sky_coef <- scan(file, "character", skip = 1)
  sky_coef <- data.frame(
    name = Map(function(x) x[1], strsplit(sky_coef, "=")) %>% unlist(),
    value = Map(function(x) x[2], strsplit(sky_coef, "=")) %>% unlist()
  )
  sky_coef[c(2, 1, 5, 4, 3), 2] %>% sub(",", ".", .) %>% as.numeric()
}
