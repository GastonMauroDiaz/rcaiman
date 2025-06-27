#' compute complementary gradients
#'
#' @inheritParams expand_noncircular
#'
#' @return An object of class [SpatRaster-class]
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' com <- compute_complementary_gradients(caim)
#' mem <- max(com$blue_yellow, com$cyan_red)
#'
#' thrs <- calc_thrs(mem, z, a,  m,
#'                   angle_width = 30,
#'                   fov = 60,
#'                   method = "thr_twocorner_uc")
#'
#' thrs_f <- filter_thrs(thrs, r, z, a, 2.5, FALSE, n_min = 0)
#'
#' plot(apply_thr(mem, mean(thrs[,3])))
#' }
compute_complementary_gradients <- function(caim) {
  stopifnot(class(caim) == "SpatRaster")
  stopifnot(all(names(caim) == c("Red", "Green", "Blue")))

  # caim <- terra::focal(caim, 3, mean)

  R <- caim$Red
  G <- caim$Green
  B <- caim$Blue
  brightness <- R + G + B

  thr <- stats::quantile(brightness[brightness != 0], 0.1, na.rm = TRUE)
  low_exposure <- terra::rast(R)
  low_exposure[] <- stats::plogis(brightness[],
                                  thr,
                                  stats::IQR(brightness[], na.rm = TRUE))

  magenta_green   <- (R - G + B) / brightness * low_exposure
  blue_yellow <- (-R - G + B) / brightness * low_exposure
  cyan_red  <- (-R + G + B) / brightness * low_exposure

  names(magenta_green) <- "magenta_green"
  names(blue_yellow) <- "blue_yellow"
  names(cyan_red) <- "cyan_red"
  c(magenta_green, blue_yellow, cyan_red)
}

