#' Local fuzzy thresholding
#'
#' This function was first presented in \insertCite{Diaz2015;textual}{rcaiman}.
#' It uses a threshold value as the location parameter of a logistic membership
#' function whose scale parameter depends on a variable, here named `mem`. This
#' dependence can be explained as follows: if the variable is equal to `1`, then
#' the membership function is same as a threshold function because the scale
#' parameter is `0`; lowering the variable increases the scale parameter, thus
#' blurring the threshold because it decreases the steepness of the curve. Since
#' the variable is defined pixel by pixel, this should be considered as a
#' **local** fuzzy thresholding method.
#'
#' Argument `m` can be used to affect the automatic estimation of `thr` and
#' `fuzziness`.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`).
#'
#'
#' @param lightness [SpatRaster-class]. A normalized greyscale image (see
#'   [normalize_minmax()]).
#' @inheritParams masking
#' @param mem [SpatRaster-class]. It is the scale parameter of the logistic
#'   membership function. Typically it is obtained with [membership_to_color()].
#' @param thr Numeric vector of length one. Location parameter of the logistic
#'   membership function. Use `NULL` to estimate it automatically with
#'   [thr_isodata()].
#' @param fuzziness Numeric vector of length one. This number is a constant
#'   value that scales `mem`. Use `NULL` to estimate it automatically as the
#'   midpoint between the maximum and minimum values of `lightness`.
#'
#' @references \insertAllCited{}
#'
#' @return An object of class [SpatRaster-class] with same pixel dimensions than
#'   `caim`. Depending on `mem`, changes could be subtle.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(2132/2, lens("Olloclip"))
#' a <- azimuth_image(z, orientation = 180)
#' zenith_colrow <- c(1063, 771)/2
#' caim <- expand_noncircular(caim, z, zenith_colrow)
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#' caim <- normalize_minmax(caim)
#'
#' mem <- membership_to_color(caim, sRGB(0.1, 0.4, 0.8))
#' plot(mem)
#' mem_thr <- local_fuzzy_thresholding(caim$Blue, m,
#'                                     mem$membership_to_grey,
#'                                     fuzziness = 0.1)
#' plot(mem_thr)
#' }
local_fuzzy_thresholding <- function (lightness,
                                      m,
                                      mem,
                                      thr = NULL,
                                      fuzziness = NULL) {
  .is_single_layer_raster(lightness, "lightness")
  .was_normalized(lightness, "lightness")
  .is_single_layer_raster(m, "m")
  .is_logic_and_NA_free(m, "m")
  .is_single_layer_raster(mem, "mem")
  terra::compareGeom(lightness, m)
  terra::compareGeom(lightness, mem)

  if (is.null(thr)) {
    thr <- thr_isodata(lightness[m] %>% as.numeric())
  }
  if (is.null(fuzziness)) {
    fuzziness <- (max(lightness[m]) - min(lightness[m])) / 2
  }
  fun <- function(lightness, mem) {
    suppressWarnings(stats::plogis(lightness, thr, fuzziness * mem))
  }
  mem <- fun(terra::values(lightness), terra::values(mem))
  terra::values(lightness) <- mem
  names(lightness) <- "membership_to_values_above_the_threshold_value"
  lightness
}
