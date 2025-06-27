#' Title
#'
#' @inheritParams apply_thr
#' @inheritParams sky_grid_segmentation
#' @inheritParams calc_co
#' @param method Character vector of length one.
#'
#' @returns data.frame
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
#' thrs <- calc_thrs(r, z, a,  m,
#'                   angle_width = 30,
#'                   fov = 60,
#'                   method = "thr_isodata")
#'
#' thrs_f <- filter_thrs(thrs, r, z, a, 2.5, TRUE, n_min = 3)
#' thr <- interpolate_planar(thrs_f$rr$sky_points, r,
#'                           rmax = ncol(r)/7, col_id = "dn")
#' plot(thr)
#' plot(apply_thr(r, thr))
#' }
calc_thrs <- function(r, z, a, m,
                      angle_width = 6,
                      fov = 30,
                      method = "thr_twocorner_uc") {

  sky_points <- sky_grid_points(z, a, angle_width)
  rr <- extract_rel_radiance(r, z, a, sky_points, NULL, FALSE)

  .fun <- switch (method,
    thr_isodata = thr_isodata,
    thr_twocorner_tl = function(x) thr_twocorner(x, window_length = rep(1,11,2))$tl,
    thr_twocorner_tm = function(x) thr_twocorner(x, window_length = rep(1,11,2))$tm,
    thr_twocorner_th = function(x) thr_twocorner(x, window_length = rep(1,11,2))$th,
    thr_twocorner_lc = function(x) thr_twocorner(x, window_length = rep(1,11,2))$lc,
    thr_twocorner_uc = function(x) thr_twocorner(x, window_length = rep(1,11,2))$uc
  )

  sky_points2 <- sky_grid_points(z, a, angle_width/2^2)
  thrs <- NULL
  rthr <- rast(r)
  for (i in 1:nrow(sky_points)) {
    rthr[] <- NA
    m_fov <- select_circumsolar_region(z, a,
                                       rr$sky_points[i, c("z", "a")] %>%
                                         as.numeric(),
                                       fov / 2)
    m_fov[!m] <- 0
    tryCatch(rthr[m_fov] <- .fun(r[m_fov]),
             error = function(e) NA)
    thr <- extract_dn(rthr, sky_points2, use_window = FALSE)
    thr <- thr[!is.na(thr[,3]),]
    thrs <- rbind(thrs, thr)
  }

  names(thrs)[3] <- method
  thrs
}
