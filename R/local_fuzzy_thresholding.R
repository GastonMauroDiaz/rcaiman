#' local fuzzy thresholding
#'
#' This function was first presented in \insertCite{Diaz2015;textual}{rcaiman}.
#' It uses a threshold value as the location parameter of a logistic membership
#' function whose scale parameter depends on a variable, here named \code{mem}.
#' This dependence can be explained as follows: if the variable is equal to
#' \code{1}, then the membership function is same as a threshold function
#' because the scale parameter is \code{0}; lowering the variable increases the
#' scale parameter, thus blurring the threshold because it decreases the
#' steepness of the curve. Since the variable is defined pixel by pixel, this
#' should be considered as a \strong{local} fuzzy thresholding method.
#'
#' Argument \code{m} can be used to affect the automatic estimation of
#' \code{thr} and \code{fuzziness}.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package.
#'
#'
#' @param lightness \linkS4class{SpatRaster}. A normalized greyscale image (see
#'   \code{\link{normalize}}).
#' @param m \linkS4class{SpatRaster}. A mask. For hemispherical photographs,
#'   check \code{\link{mask_hs}}.
#' @param mem \linkS4class{SpatRaster}. It is the scale parameter of the
#'   logistic membership function. Typically it is obtained with
#'   \code{\link{membership_to_color}}.
#' @param thr Numeric vector of length one. Location parameter of the logistic
#'   membership function. Use \code{NULL} to estimate it automatically with
#'   \code{\link{thr_isodata}}.
#' @param fuzziness Numeric vector of length one. This number is a constant
#'   value that scales \code{mem}. Use \code{NULL} to estimate it automatically
#'   as the midpoint between the maximum and minimum values of \code{lightness}.
#'
#' @references \insertAllCited{}
#'
#' @return An object of class \linkS4class{SpatRaster} with same pixel
#'   dimensions than \code{caim}. Depending on \code{mem}, changes could be
#'   subtle; however, they should be in the direction of showing more contrast
#'   between the sky and plant pixels than any of the individual bands from
#'   \code{caim}.
#'
#'
#' @export
#' @family Pre-processing Functions
#' @examples
#' \donttest{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' target_color <- sRGB(matrix(c(0.529, 0.808, 0.921), ncol = 3))
#' mem <- membership_to_color(caim, target_color)
#' m <- !is.na(z)
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
