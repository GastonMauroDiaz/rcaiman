#' Read manual input
#'
#' Read manual input stored in an HSP project.
#'
#' Refer to the Details section of function
#' [write_sky_points()].
#'
#' @inheritParams write_sky_points
#' @family HSP Functions
#'
#' @return A list of numeric vectors named *weight*, *max_points*,
#'   *angle*, *point_radius*, *sun_coord*, *sky_points* and
#'   *zenith_dn.*
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
       sun_coord = sun_mark,
       sky_points = sky_marks,
       zenith_dn = zenith_dn)
}
