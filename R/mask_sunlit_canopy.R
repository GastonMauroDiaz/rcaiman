#' Mask sunlit canopy
#'
#' It is a wrapper function around [membership_to_color()]. It was developed
#' with images in sRGB color space \insertCite{Diaz2023}{rcaiman}.
#'
#' @inheritParams enhance_caim
#'
#' @return An object of class [SpatRaster-class] with values `0` and `1`.
#'
#' @export
#'
#' @family Segmentation Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' plotRGB(caim)
#' caim <- normalize(caim)
#' m <- mask_sunlit_canopy(caim)
#' plot(m)
#' }
mask_sunlit_canopy <- function(caim, m = NULL) {
  if (is.null(m)) m <- !is.na(caim$Red)
  mem <- membership_to_color(caim, sRGB(0.25,0.75,0))
  sunlit_canopy <- mem$membership_to_target_color > 0.75
  sunlit_canopy_pct <- (sum(sunlit_canopy[], na.rm = TRUE) /
                          sum(m[], na.rm = TRUE)) * 100
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
