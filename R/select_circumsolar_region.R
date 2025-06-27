#' Select circumsolar region
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams fit_cie_sky_model
#' @inheritParams extract_sun_zenith_azimuth
#'
#' @return An object of class [SpatRaster-class] with values `0` and
#'   `1`.
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
#' sun_zenith_azimuth <- extract_sun_zenith_azimuth(r, z, a, m, NULL)
#' m_sun <- select_circumsolar_region(z, a, sun_zenith_azimuth, 20)
#' plot(m_sun)
#' }
select_circumsolar_region <- function(z, a, sun_zenith_azimuth, chi_max_sun) {
  m_sun <- rast(z)
  m_sun[] <- calc_spherical_distance(z[] %>% .degree2radian(),
                                     a[] %>% .degree2radian(),
                                     sun_zenith_azimuth[1] %>% .degree2radian(),
                                     sun_zenith_azimuth[2] %>% .degree2radian())
  m_sun <- !apply_thr(m_sun, .degree2radian(chi_max_sun))
  m_sun[is.na(z)] <- 0
  m_sun
}
