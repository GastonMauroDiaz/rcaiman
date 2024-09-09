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
#'   [normalize()]).
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
#' @family Pre-processing Functions
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' caim <- normalize(caim)
#'
#' # ImageJ can be used to digitize points
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' img_points <- read.csv(path)
#' img_points <- img_points[c("Y", "X")]
#' colnames(img_points) <- c("row", "col")
#' head(img_points)
#' target_color <- extract_dn(caim, img_points, fun = median)
#' as(target_color, "HSV")
#' target_color <- HSV(240, 0.85, 0.5) #to increase saturation
#'
#' mem <- membership_to_color(caim, target_color)
#' mem_thr <- local_fuzzy_thresholding(mean(caim), m,  mem$membership_to_grey)
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
  names(lightness) <- "membership_to_values_above_the_threshold"
  lightness
}
