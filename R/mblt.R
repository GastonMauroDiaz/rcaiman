#' Threshold image
#'
#' Transform background digital number into threshold values.
#'
#' Transform background digital number into threshold values by means of the
#' Equation 1 presented in \insertCite{Diaz2018}{rcaiman}, which is a linear
#' function with the slope modified by a weighting parameter. This simple
#' function was developed by studying canopy models, also known as targets,
#' which are planes with holes made of a rigid and dark material. These models
#' were backlighted with homogeneous lighting and photographs were acquired with
#' a Nikon Coolpix 5700 in JPEG format and then gamma back corrected with a
#' default gamma value equal to 2.2 (see \code{\link{gbc}}). Results clearly
#' shown that the optimal threshold value was linearly related with the
#' background digital number. So, encouraging to shift attention from finding
#' the optimal threshold value to obtain the digital number (DN) of the
#' background, which is usually the sky. Functions \code{\link{model_sky_dn}}
#' and \code{\link{fit_trend_surface_to_sky_dn}} address that topic. In other
#' words, \code{thr_image} function should be used in combination with those
#' functions. An out-of-the-box solution based on these principles is available
#' through the \code{\link{mblt}} function.
#'
#' It is worth noting that Equation 1 was developed with 8-bit images, so
#' calibration of new coefficient should be done in the 0 to 255 domain since
#' that is what \code{\link{thr_image}} expect, although the input \code{dn}
#' should be normalized. The latter was a design decision aiming to harmonize
#' the whole package.
#'
#' @param dn Numeric vector or \linkS4class{RasterLayer}. Digital number of the
#'   background. These values should be normalized and, if they are extracted
#'   from JPEG image, gamma back corrected.
#' @param w Numeric vector of length one. Weight affecting the slope
#'   coefficient. See \insertCite{Diaz2018}{rcaiman}.
#' @param type Character vector of length one. Default is "Generic". Currently,
#'   the only available calibrated values are from
#'   \insertCite{Diaz2018}{rcaiman}. Use "Nikon_Coolpix_5700" to use them.
#' @param intercept,slope Numeric vector of length one. Default is NULL. These
#'   arguments allow providing curstomized coefficients. See
#'   \insertCite{Diaz2018}{rcaiman} for details.
#'
#' @export
#'
#' @family mblt functions
#' @seealso \code{\link{normalize}}, \code{\link{gbc}}, \code{\link{apply_thr}}
#'   and \code{\link{regional_thresholding}}.
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' thr_image(125)
thr_image <- function (dn,
                       w = 0.5,
                       type = "Generic",
                       intercept = NULL,
                       slope = NULL) {
  stopifnot(length(w) == 1)
  stopifnot(class(w) == "numeric")
  stopifnot(class(type) == "character")
  if (!is.null(intercept)) stopifnot(length(intercept) == 1)
  if (!is.null(slope)) stopifnot(length(slope) == 1)

  dn <- dn * 255

  if (all(is.null(intercept), is.null(slope))) {

    my_coef <- switch(type,
                      Generic = c(-8, 1),
                      Nikon_Coolpix_5700 = c(-7.7876, 0.9485))

    intercept <- my_coef[1]
    slope <- my_coef[2]
  }

  if (length(dn[dn > 255]) > 0) {dn[dn > 255] <- 255}
  thr <- intercept  + slope * dn * w
  if (length(thr[thr < 0]) > 0) {thr[thr < 0] <- 0}
  thr / 255
}



#' Model-based local thresholding
#'
#' Out-of-the-box version of the model-based local thresholding (MBLT)
#' algorithm.
#'
#' This function is a hard-coded version of a MBLT pipeline that starts with a
#' working binarized image and ended with a refined binarized image. The
#' pipeline combines \code{\link{regional_thresholding}},
#' \code{\link{model_sky_dn}}, \code{\link{fit_trend_surface_to_sky_dn}}, and
#' \code{\link{thr_image}}. The code can be easily inspected by calling
#' \code{mblt}. Advanced users could use that code as a template.
#'
#' The MBLT algorithm was first presented in \insertCite{Diaz2018}{rcaiman}.
#' This version differs from that in the following main aspects:
#'
#' \itemize{
#'
#' \item This version implements a regional threholding approach as first step
#' instead of a global one.
#'
#' \item It does not use asynchronous acquisition under the open sky.
#'
#' \item The working binarized image is not refined in the second step by using
#' the sky digital number (sDN) obtained with Equation 4
#' \code{\link{model_sky_dn}} and the local threshold values calculated with
#' Equation 1 \code{\link{thr_image}}. Instead, the sDN are used as filling
#' source for trend surface fitting \code{\link{fit_trend_surface_to_sky_dn}}.
#'
#' \item The sDN obtained by trend surface fitting is merged with the sDN
#' obtained with Equation 4. To merge them, a weighted average is calculated, in
#' which weight are proportional to the zenith angle. Near the zenith the values
#' obtained by means of trend surface fitting prevail over the ones obtained
#' with Equation 7, and the opposite occur near the horizon.
#'
#' \item This version use a generic Equation 1 \code{\link{thr_image}} instead
#' of a specific one.
#'
#' }
#'
#'
#' @inheritParams model_sky_dn
#' @inheritParams thr_image
#'
#' @return If \code{w} has length equal to one, \linkS4class{RasterLayer}. Else,
#'   \linkS4class{RasterStack}.
#'
#' @export
#' @family mblt functions
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' \dontrun{
#' path <- getwd()
#' my_file <- paste0(path, "/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' bin <- mblt(blue, z, a)
#' }
mblt <- function(r, z, a, w = 0.5) {
  seg <- sky_grid_segmentation(z, a, 30)
  m <- mask_image(z, zlim = c(0,70))

  prob <- 1
  sky_m <- NA
  while (length(sky_m) == 1) {
    if (prob < 0.9) {
      stop(paste("The function is not working properly.",
                 "The problem might be related to inputs.",
                 "Please, make sure they are OK."))
    }
    prob <- prob - 0.01
    bin <- regional_thresholding(r, seg, "Diaz2018", 0.9, "Generic", prob)
    sky_m <- model_sky_dn(r, z, a, bin)
  }
  sky_m <- sky_m$image

  sky_s <- fit_trend_surface_to_sky_dn(r, z, m, bin, sky_m)$image

  mask <- (sky_m - sky_s) > 0.5
  sky_s[mask] <- NA

  sky_combo <- sky_s * z / 90 + sky_m * (1 - z / 90)

  sky <- cover(sky_combo, sky_m)

  fun <- function(w) {
    apply_thr(r, thr_image(sky, w))
  }

  if (is.null(w)) {
    return(sky_s)
  } else {
    if (length(w) == 1) {
      bin <- fun(w)
    } else {
      bin <- Map(fun, w)
      bin <- stack(bin)
      names(bin) <- paste0("w", w)
    }
    return(bin)
  }
}






