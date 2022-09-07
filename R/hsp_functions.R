#' Write sky points
#'
#' Create a special file to interface with the HSP software.
#'
#' This function is part of a workflow that connects this package with the HSP
#' software package \insertCite{Lang2013}{rcaiman}.
#'
#' This function was designed to be called after
#' \code{\link{extract_sky_points}}. The \code{r} argument provided to
#' \code{\link{extract_sky_points}} should be an image pre-processed by the HSP
#' software. Those images are stored as PGM files in the subfolder "manipulate"
#' of the project folder (which will be in turn a subfolder of the
#' "project\strong{s}" folder). Those PGM files can be read with
#' \code{\link{read_caim}}.
#'
#' The \code{img_name} argument of \code{write_sky_points()} should be the name
#' of the file associated to the aforementioned \code{r} argument.
#'
#' The following code exemplifies how this package can be used in conjunction
#' with HSP software.  The code assumes that the user is working within an
#' RStudio project located in the HSP project folder.
#'
#' \preformatted{
#' r <- read_caim("manipulate/IMG_1014.pgm") %>% normalize()
#' plot(r)
#' z <- zenith_image(ncol(r), lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' mblt <- ootb_mblt(r, z, a)
#' bin <- find_sky_pixels_nonnull(r, mblt$sky_s, g)
#' bin <- mask_hs(z, 0, 85) & bin
#'
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' write_sun_coord(sun_coord$row_col, ".", "IMG_1014")
#'
#' sky_points <- extract_sky_points(r, bin, g)
#' write_sky_points(sky_points, ".", "IMG_1014")
#' }
#'
#' @param sky_points An object of the class \emph{data.frame}. The result of a
#'   calling to \code{\link{extract_sky_points}}.
#' @param path_to_HSP_project Character vector of length one. Path to the HSP
#'   project folder. For instance,
#'   "C:/Users/johndoe/Documents/HSP/Projects/my_prj/".
#' @param img_name Character vector of length one. For instance, "DSCN6342.pgm"
#'   or "DSCN6342". See details.
#'
#' @family HSP Functions
#'
#' @references \insertAllCited{}
#'
#' @return None. A file will be written in the HSP project folder.
#'
#' @export
write_sky_points <- function(sky_points, path_to_HSP_project, img_name) {
  no <- nrow(sky_points)

  col.row_coordinates <- paste(sky_points$col, sky_points$row, "3", sep = ".")
  col.row_coordinates <- paste(col.row_coordinates, collapse = " ")

  sky_points <- c(no, col.row_coordinates)
  sky_points <- data.frame(sky_points)

  img_name <- filenamer::as.filename(img_name)
  img_name <- filenamer::trim_ext(img_name) %>% as.character()

  utils::write.table(sky_points, file.path(path_to_HSP_project,
                                           "manipulate",
                                           paste0(img_name, "_points.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}


#' Write sun coordinates
#'
#' Create a special file to interface with the HSP software.
#'
#' Refer to the Details section of function \code{\link{write_sky_points}}.
#'
#' @param sun_coord Numeric vector of length two. Raster coordinates of the
#'   solar disk that can be obtained by calling to
#'   \code{\link{extract_sun_coord}}. \strong{TIP}: if the output of
#'   \code{extrac_sun_coord()} is \code{sun_coord}, then you should provide here
#'   this: \code{sun_coord$row_col}. See also
#'   \code{\link{row_col_from_zenith_azimuth}} in case you want to provide
#'   values based on date and time of acquisition and the R package 'suncalc'.
#' @inheritParams write_sky_points
#'
#'
#' @family HSP Functions
#'
#' @return None. A file will be written in the HSP project folder.
#'
#' @export
write_sun_coord <- function(sun_coord, path_to_HSP_project, img_name) {
  sun_coord <- paste(sun_coord[c(2,1)], collapse = ".")

  img_name <- filenamer::as.filename(img_name)
  img_name <- filenamer::trim_ext(img_name) %>% as.character()

  utils::write.table(sun_coord, file.path(path_to_HSP_project,
                                         "manipulate",
                                         paste0(img_name, "_sun.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}

#' Read manual input
#'
#' Read manual input stored in an HSP project.
#'
#' Refer to the Details section of function
#' \code{\link{write_sky_points}}.
#'
#' @inheritParams write_sky_points
#' @family HSP Functions
#'
#' @return A list of numeric vectors named \emph{weight}, \emph{max_points},
#'   \emph{angle}, \emph{point_radius}, \emph{sun_mark}, \emph{sky_marks} and
#'   \emph{zenith_dn.}
#'
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
#' Refer to the Details section of function
#' \code{\link{write_sky_points}}.
#'
#' @inheritParams write_sky_points
#' @family HSP Functions
#' @return Numeric vector of length five.
#' @seealso \code{\link{cie_sky_model_raster}}
#'
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

#' Row and col numbers from zenith and azimuth angles
#'
#' @inheritParams ootb_mblt
#' @inheritParams zenith_image
#' @param za Numeric vector of length two. Zenith and azimuth angles in degrees.
#'
#' @export
#'
#' @family HSP Functions
#' @return Numeric vector of length two.
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' row_col_from_zenith_azimuth(z, c(45, 270), lens())
row_col_from_zenith_azimuth <- function(r, za, lens_coef) {
  .is_single_layer_raster(r, "r")
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
  p <- terra::vect(matrix(cart[1:2], ncol = 2))
  terra::crs(p) <- terra::crs(r)
  e <- terra::ext(r)
  terra::ext(r) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
  ir <- terra::rasterize(p, r)
  l <- list(za, terra::cells(r, terra::ext(p)) %>%
                   terra::rowColFromCell(r, .) %>% as.numeric())
  names(l) <- c("zenith_azimuth", "row_col")
  l
}


#' Zenith and azimuth angles from row and col numbers
#'
#' @inheritParams ootb_mblt
#' @param row_col Numeric vector of length two. Row and col numbers.
#' @inheritParams zenith_image
#'
#' @export
#'
#' @family HSP Functions
#' @examples
#' z <- zenith_image(1000, lens_coef = lens())
#' zenith_azimuth_from_row_col(z, c(501, 750), lens())
zenith_azimuth_from_row_col <- function(r, row_col, lens_coef) {
  .is_single_layer_raster(r, "r")
  stopifnot(ncol(r) == nrow(r))
  stopifnot(is.numeric(lens_coef))
  stopifnot(is.numeric(row_col))
  stopifnot(length(row_col) == 2)

  #get azimuth
  e <- terra::ext(r)
  terra::ext(r) <- terra::ext(-pi/2,pi/2,-pi/2,pi/2)
  xy <- terra::cellFromRowCol(r, row_col[1], row_col[2]) %>%
    terra::xyFromCell(r, .)
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
  l <- list(c(zenith, azimuth), row_col)
  names(l) <- c("zenith_azimuth", "row_col")
  l
}


