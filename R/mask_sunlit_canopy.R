#' Mask sunlit canopy
#'
#' @inheritParams enhance_caim
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1} that masks pixels that are very likely sunlit canopy.
#'
#' @export
#'
#' @examples
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' m<- is.na(z)
#' sunlit_canopy <- mask_sunlit_canopy(caim, m)
#' plot(sunlit_canopy)
mask_sunlit_canopy <- function(caim, m) {
  mem <- membership_to_color(caim, sRGB(0.25,0.75,0))
  sunlit_canopy <- mem$membership_to_target_color > 0.75
  sunlit_canopy_pct <- (sum(sunlit_canopy[]) / sum(m[])) * 100
  if (sunlit_canopy_pct > 0.1) {
    sunlit_canopy <- caim[sunlit_canopy] %>%
      apply(., 2, median) %>%
      as.numeric() %>%
      matrix(., ncol = 3)
    mem <- membership_to_color(caim, sRGB(sunlit_canopy))
    sunlit_canopy <- mem$membership_to_target_color > 0.75
  }
  sunlit_canopy
}
