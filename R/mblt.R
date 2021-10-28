#' Threshold image
#'
#' Transform background digital number into threshold values.
#'
#' Transform background digital number into threshold values by means of the
#' Equation 1 presented in \insertCite{Diaz2018;textual}{rcaiman}, which is a
#' linear function with the slope modified by a weighting parameter. This simple
#' function was found by studying canopy models, also known as targets, which
#' are planes with holes made of a rigid and dark material. These models were
#' backlighted with homogeneous lighting, photographed with a Nikon Coolpix 5700
#' set to acquire in JPEG format, and those images were gamma back corrected
#' with a default gamma value equal to 2.2 (see \code{\link{gbc}}). Results
#' clearly shown that the optimal threshold value was linearly related with the
#' background digital number. Therefore, that shifts the aim from finding the
#' optimal threshold to obtaining the background DN as if the canopy were not
#' there. Functions \code{\link{fit_cone_shaped_model}} and
#' \code{\link{fit_trend_surface_to_sky_dn}} address that topic.
#'
#' It is worth noting that Equation 1 was developed with 8-bit images, so
#' calibration of new coefficient should be done in the 0 to 255 domain since
#' that is what \code{\link{thr_image}} expect, although the input \code{dn}
#' should be normalized. The latter --that might sound counter intuitive-- was a
#' design decision aiming to harmonize the whole package.
#'
#' To apply the weighting parameter (w) from Equation 1, just provide the
#' argument \code{slope} as \code{slope_value * w}.
#'
#' Type \code{thr_image} --no parenthesis-- in the console to inspect the code,
#' which is very simple.
#'
#' @param dn Numeric vector or \linkS4class{RasterLayer}. Digital number of the
#'   background. These values should be normalized and, if they are extracted
#'   from JPEG image, gamma back corrected.
#' @param intercept,slope Numeric vector of length one. These are linear
#'   function coefficients. Please, see the Details section of
#'   \code{\link{thr_image}}.
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
#' thr_image(gbc(125), -8, 1)
thr_image <- function (dn, intercept, slope) {
  stopifnot(length(intercept) == 1)
  stopifnot(class(intercept) == "numeric")
  stopifnot(length(slope) == 1)
  stopifnot(class(slope) == "numeric")

  dn <- dn * 255

  if (.get_max(dn) > 255) warning("\"dn\" values should be normalized")
  dn[dn > 255] <- 255

  thr <- intercept  + slope * dn
  thr[thr < 0] <- 0

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
#' \code{\link{fit_cone_shaped_model}}, \code{\link{fit_trend_surface_to_sky_dn}}, and
#' \code{\link{thr_image}}. The code can be easily inspected by calling
#' \code{mblt}. Advanced users could use that code as a template.
#'
#' The MBLT algorithm was first presented in
#' \insertCite{Diaz2018;textual}{rcaiman}. This version differs from that in the
#' following main aspects:
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
#' \code{\link{fit_cone_shaped_model}} and the local threshold values calculated with
#' Equation 1 \code{\link{thr_image}}. Instead, the sDN are used as filling
#' source for trend surface fitting \code{\link{fit_trend_surface_to_sky_dn}}.
#'
#' \item The sDN obtained by trend surface fitting is merged with the sDN
#' obtained with Equation 4. To merge them, a weighted average is calculated,
#' being weights proportional to the zenith angle. Near the zenith, values
#' obtained by means of trend surface fitting prevail over the ones obtained
#' with Equation 7, and the opposite occur near the horizon.
#'
#' }
#'
#'
#' @inheritParams fit_cone_shaped_model
#' @inheritParams thr_image
#' @param w Numeric vector. Weighting parameter for \code{slope}.
#'
#'
#' @return If \code{w} has length equal to one, \linkS4class{RasterLayer}. Else,
#'   \linkS4class{RasterStack}. If \code{w} is \code{NULL}, the above canopy
#'   reconstruction is returned.
#'
#' @export
#' @family mblt functions
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
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
mblt <- function(r, z, a, intercept = -8, slope = 1, w = 0.5) {
  .check_if_r_was_normalized(r)
  seg <- sky_grid_segmentation(z, a, 30)

  prob <- 1
  sky_m <- NA
  while (length(sky_m) == 1) {
    if (prob < 0.9) {
      stop(paste("The function is not working properly.",
                 "The problem might be related to inputs.",
                 "Please, make sure they are OK."))
    }
    prob <- prob - 0.01
    bin <- regional_thresholding(r, seg, "Diaz2018", 0, 1, prob)
    sky_m <- fit_cone_shaped_model(r, z, a, bin)
  }
  sky_m <- sky_m$image

  m <- mask_image(z, zlim = c(0,70))
  sky_s <- fit_trend_surface_to_sky_dn(r, z, m, bin, sky_m)$image

  mask <- (sky_m - sky_s) > 0.5
  sky_s[mask] <- NA

  sky_combo <- sky_s * z / 90 + sky_m * (1 - z / 90)

  sky <- cover(sky_combo, sky_m)

  if (is.null(w)) {
    return(sky)
  } else {
    if (length(w) == 1) {
      bin <- fun(w)
    } else {
      .bin_fun <- function(w) {
        apply_thr(r, thr_image(sky,intercept, slope * w))
      }
      bin <- Map(.bin_fun, w)
      bin <- stack(bin)
      names(bin) <- paste0("w", w)
    }
    return(bin)
  }
}






