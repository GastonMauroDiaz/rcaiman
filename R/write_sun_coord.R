#' Write sun coordinates
#'
#' Create a special file to interface with the HSP software.
#'
#' Refer to the Details section of function [write_sky_points()].
#'
#' @param sun_row_col Numeric vector of length two. Raster coordinates (row and
#'   column) of the solar disk. See
#'   also [row_col_from_zenith_azimuth()] in case you want to provide values
#'   based on date and time of acquisition and the `suncalc` package.
#' @inheritParams write_sky_points
#'
#' @seealso [row_col_from_zenith_azimuth()]
#'
#' @return None. A file will be written in the HSP project folder.
#'
#' @export
write_sun_coord <- function(sun_row_col, path_to_HSP_project, img_name) {
  sun_col_row <- paste(sun_row_col[c(2,1)], collapse = ".")

  img_name <- filenamer::as.filename(img_name)
  img_name <- filenamer::trim_ext(img_name) %>% as.character()

  utils::write.table(sun_col_row, file.path(path_to_HSP_project,
                                            "manipulate",
                                            paste0(img_name, "_sun.conf")),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              fileEncoding = "UTF-8", eol = "\n")
}
