#' local fuzzy thresholding
#'
#' This function is presented in \insertCite{Diaz2015;textual}{rcaiman}. It uses
#' a threshold value as the location parameter of a logistic membership function
#' whose scale parameter depends on a variable, here named \code{mem}. This
#' dependence can be explained as follows: if the variable is equal to \code{1},
#' then the membership function is same as a threshold function because the
#' scale parameter is \code{0}; lowering the variable increases the scale
#' parameter, thus blurring the threshold because it decreases the steepness of
#' the curve. Since the variable is defined pixel by pixel, this should be
#' considered as a \strong{local} fuzzy thresholding method.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015}{rcaiman}.
#'
#'
#' @param lightness \linkS4class{RasterLayer}. A normalized greyscale image, the
#'   lightness value. Values should range between zero and one --please see
#'   \code{\link{normalize}}.
#' @inheritParams fit_trend_surface
#' @param mem \linkS4class{RasterLayer}. It is the scale parameter of the
#'   logistic membership function. Typically it is obtained with
#'   \code{\link{membership_to_color}}.
#' @param thr Numeric vector of length one. Location parameter of the logistic
#'   membership function. Use \code{NULL} (default) to estimate it
#'   automatically with the function \code{\link[autothresholdr]{auto_thresh}},
#'   method \code{"IsoData"}.
#' @param fuzziness Numeric vector of length one. This number is a constant that
#'   scale \code{mem}. Use \code{NULL} (default) to estimate it
#'   automatically as the midpoint between the maximum and minimum values of
#'   \code{lightness}.
#'
#' @references \insertAllCited{}
#'
#' @details Argument \code{m} can be used to affect the estimation of \code{thr}
#'   and \code{fuzziness}.
#'
#' @export
#' @family Fuzzy logic functions
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' target_color <- colorspace::sRGB(matrix(c(0.529, 0.808, 0.921), ncol = 3))
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
  .check_if_r_was_normalized(lightness, "lightness")
  if (!compareRaster(lightness, m, stopiffalse = FALSE)) {
    stop("\"x\" should match pixel by pixel whit \"m\".")
  }

  if (is.null(thr)) {
    if (!requireNamespace("autothresholdr", quietly = TRUE)) {
      stop(paste(
        "Package \"autothresholdr\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
      )
    }
    dns <- lightness[m]
    thr <- autothresholdr::auto_thresh(round(dns * 255), "IsoData")[1] / 255
  }
  if (is.null(fuzziness)) {
    fuzziness <- (max(lightness[m]) - min(lightness[m])) / 2
  }

  mem <- overlay(lightness, mem, fun = function(lightness, mem) {
    suppressWarnings(stats::plogis(lightness, thr, fuzziness * (1 - mem)))
  })
  names(mem) <- "membership_to_values_above_the_threshold"
  mem
}
