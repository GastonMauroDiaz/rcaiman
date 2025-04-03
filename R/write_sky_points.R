#' Write sky points
#'
#' Create a special file to interface with the HSP software.
#'
#' This function is part of a workflow that connects this package with the HSP
#' software package \insertCite{Lang2013}{rcaiman}.
#'
#' This function was designed to be called after
#' [extract_sky_points()]. The `r` argument provided to
#' [extract_sky_points()] should be an image pre-processed by the HSP
#' software. Those images are stored as PGM files in the subfolder "manipulate"
#' of the project folder (which will be in turn a subfolder of the
#' "project**s**" folder). Those PGM files can be read with
#' [read_caim()].
#'
#' The `img_name` argument of `write_sky_points()` should be the name
#' of the file associated to the aforementioned `r` argument.
#'
#' The following code exemplifies how this package can be used in conjunction
#' with the HSP software. The code assumes that the user is working within an
#' RStudio project located in the HSP project folder.
#'
#' \preformatted{
#' r <- read_caim("manipulate/IMG_1014.pgm") %>% normalize()
#' plot(r)
#' z <- zenith_image(ncol(r), lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' mblt <- ootb_mblt(r, z, a)$bin
#' bin <- select_sky_vault_region(z, 0, 85) & bin
#'
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' write_sun_coord(sun_coord$row_col, ".", "IMG_1014")
#'
#' sky_points <- extract_sky_points(r, bin, g)
#' write_sky_points(sky_points, ".", "IMG_1014")
#' }
#'
#' @param sky_points An object of the class *data.frame*. The result of a
#'   calling to [extract_sky_points()].
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
