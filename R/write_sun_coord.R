#' Write sun coordinates
#'
#' Create a special file to interface with the HSP software.
#'
#' Refer to the Details section of function [write_sky_points()].
#'
#' @param sun_coord Numeric vector of length two. Raster coordinates of the
#'   solar disk that can be obtained by calling to
#'   [extract_sun_coord()]. **TIP**: if the output of
#'   `extrac_sun_coord()` is `sun_coord`, then you should provide here
#'   this: `sun_coord$row_col`. See also
#'   [row_col_from_zenith_azimuth()] in case you want to provide
#'   values based on date and time of acquisition and the `suncalc` package.
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
