#' Threshold function
#'
#' Threshold value in function of the background digital number.
#'
#' TODO: Add details about the calibration range in 8bit but input normalized
#'
#' @param dn Numeric vector. Digital number of the background. The background
#'   should be lighter than the objects. In canopy photography, the typical
#'   background is the sky.
#' @param w Numeric vector of length one. Weight. See reference.
#' @param type Character vector of length one. Default is "Generic". Currently,
#'   the only available calibrated values are from
#'   \insertCite{Diaz2018}{rcaiman}. Use "Nikon_Coolpix_5700" to use them.
#' @param intercept,slope Numeric vector of length one. Default is NULL. You can
#'   provide your own calibrated values. See \insertCite{Diaz2018}{rcaiman} for
#'   details.
#'
#' @export
#'
#' @references
#' \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' thr_image(125)
thr_image <- function (dn,
                       w = 0.5,
                       type = "Generic",
                       intercept = NULL,
                       slope = NULL) {
  #stopifnot(class(dn) == "numeric")
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



mblt <- function(r, z, a, w) {
  seg <- sky_grid_segmentation(z, a, 30)
  m <- mask_image(z, zlim = c(0,70))

  prob <- 1
  sky <- NA
  while (length(sky) == 1) {
    if (prob < 0.9) {
      stop(paste("The function is not working properly.",
                 "The problem might be related with the inputs.",
                 "Please, make sure they are OK."))
    }
    prob <- prob - 0.01
    bin <- regional_thresholding(r, seg, "Diaz2018", 0.9, "Generic", prob)
    sky <- model_sky_dn(r, z, a, bin)
  }

  sky_s <- fit_surface_to_sky_dn(r, z, m, bin, sky$image)
  browser()
  mask <- (sky$image - sky_s$image) > 0.5
  sky_s$image


  sky_2 <- model_sky_dn(r, z, a, bin, use_azimuth_angle = FALSE)
  sky <- stack(sky$image, sky_2$image)




  sky <- cover(sky_s$image, sky[[2]])

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
    }
    return(bin)
  }
}






