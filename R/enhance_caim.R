.gaussian2d <- function(x, y, target_a, target_b, sigma) {
  stats::dnorm(x, target_a, sigma) * stats::dnorm(y, target_b, sigma)
}

.get_gaussian_2d_parameters <- function(target_color, sigma) {
  if (class(target_color) != "LAB") target_color <- as(target_color, "LAB")
  ma <- colorspace::coords(target_color)
  target_a <- ma[, 2]
  target_b <- ma[, 3]
  if (is.null(sigma)) sigma <- sqrt(target_a^2 + target_b^2) %>% unname()
  c(target_a, target_b, sigma)
}

#' Compute membership to a color
#'
#' Compute the degree of membership to a color with two Gaussian membership
#' functions and the dimensions \emph{A} and \emph{B} from the \emph{CIE L*a*b*}
#' color space. The lightness dimension is not considered in the calculations.
#'
#' @inheritParams expand_noncircular
#' @param target_color \linkS4class{color}.
#' @param sigma Numeric vector of length one. Use \code{NULL} (default) to order
#'   to estimate it automatically as the euclidean distance between
#'   \code{target_color} and grey in the \emph{CIE L*a*b*} color space.
#'
#'
#' @return It returns an object from the class \linkS4class{RasterBrick} or
#'   \linkS4class{RasterStack} --this will depend on the input. First layer is
#'   the membership to the target color. Second layer is the membership to grey.
#'   Both membership are calculated with same \code{sigma}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSC_2881.JPG")
#' download.file("https://osf.io/x8urg/download", my_file,
#'               method = "auto", mode = "wb"
#' )
#'
#' caim <- read_caim(my_file)
#' caim <- normalize(caim, 0, 255)
#' target_color <- sRGB(matrix(c(0.529, 0.808, 0.921), ncol = 3))
#' mem <- membership_to_color(caim, target_color)
#' plot(mem)
#' }
membership_to_color <- function(caim, target_color, sigma = NULL) {
  .is_class_from_colorspace(target_color)
  .check_if_r_was_normalized(caim, "caim")
  stopifnot(names(caim) == c("Red", "Green", "Blue"))
  color <- colorspace::sRGB(values(caim))
  if (class(color) != "LAB") color <- as(color, "LAB")
  p <- .get_gaussian_2d_parameters(target_color, sigma)
  max_z <- .gaussian2d(p[1], p[2], p[1], p[2], p[3])
  x <- colorspace::coords(color)
  mem_to_color <- .gaussian2d(x[, 2], x[, 3], p[1], p[2], p[3]) / max_z

  target_color <- sRGB(matrix(c(0, 0, 0), ncol = 3))
  sigma <- p[3]
  p <- .get_gaussian_2d_parameters(target_color, sigma)
  max_z <- .gaussian2d(p[1], p[2], p[1], p[2], sigma)
  mem_to_grey <- .gaussian2d(x[, 2], x[, 3], p[1], p[2], sigma) / max_z

  r <- raster::subset(caim, 1:2)
  r$Red <- mem_to_color
  r$Green <- mem_to_grey
  names(r) <- c("membership_to_target_color", "membership_to_grey")
  r
}

#' local fuzzy thresholding
#'
#' This function uses a threshold value as the location parameter of a logistic
#' membership function whose scale parameter depends on a variable, here named
#' \code{mem}. This dependence can be explained as follows: if the variable is
#' equal to \code{1}, then the membership function is as a threshold function
#' because the scale parameter is \code{0}; lowering the variable increases the
#' scale parameter, thus blurring the threshold because it decreases the
#' steepness of the curve. Since the variable is defined pixel by pixel, this
#' should be considered as a local fuzzy thresholding method.
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
#'   membership function. Use \code{NULL} (default) to order to estimate it
#'   automatically with the function \code{\link{autoThr}}, method
#'   \code{"IsoData"}.
#' @param fuzziness Numeric vector of length one. This number is a constant that
#'   scale \code{mem}. Use \code{NULL} (default) to order to estimate it
#'   automatically as the midpoint between the maximum and minimum values of
#'   \code{lightness}.
#'
#' @details Argument \code{m} can be used to affect the estimation of \code{thr}
#'   and \code{fuzziness}.
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' caim <- normalize(caim, 0, 255)
#' target_color <- sRGB(matrix(c(0, 0, 1), ncol = 3))
#' mem <- membership_to_color(caim, target_color)
#' mem_thr <- local_fuzzy_thresholding(mean(caim), mem$membership_to_grey)
#' }
local_fuzzy_thresholding <- function (lightness,
                                      m,
                                      mem,
                                      thr = NULL,
                                      fuzziness = NULL) {
  .check_if_r_was_normalized(lightness, "lightness")
  if (!compareRaster(caim, m, stopiffalse = FALSE)) {
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
                stats::plogis(lightness, thr, fuzziness * (1 - mem))
              })
  names(mem) <- "membership_to_threshold"
  mem
}


#' Enhance canopy image
#'
#' This function uses the color perceptual attributes to enhance the contrast
#' between the sky and plants through fuzzy classification. Color has three
#' different perceptual attributes: hue, lightness, and chroma. The algorithm
#' was developed using the following premise: the color of the sky is different
#' from the color of plants. It performs the next classification rules, here
#' expressed in natural language: clear sky is blue and clouds decrease its
#' chroma; if clouds are highly dense, then the sky is achromatic, and, in such
#' cases, it can be light or dark; everything that does not match this
#' description is not sky. These linguistic rules were translated to math
#' language by means of fuzzy logic.
#'
#' @inheritParams expand_noncircular
#' @inheritParams local_fuzzy_thresholding
#' @param w_red Numeric vector of length one. Weight of the red channel. A
#'   single layer image is calculated as a weighted average of the blue and red
#'   channels. This layer is used as lightness information. The weight of the
#'   blue is the complement of the \code{w_red}.
#' @param sky_blue \linkS4class{color}. Is the \code{target_color} argument to
#'   be passed to \code{\link{membership_to_color}}.
#'
#' @details This is a pixel-wise methodology that evaluates the possibility for
#'   a pixel to be member of the class Gap. High score could mean either high
#'   membership to \code{sky_blue} or, in the case of achromatic pixels, a high
#'   membership to \code{thr}. The algorithm internally uses
#'   \code{\link{membership_to_color}} and \code{\link{fuzzy_lightness}}. The
#'   argument \code{sky_blue} is the \code{target_color} of the former function,
#'   which output is the argument \code{mem} of the latter function.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' m <- !is.na(z)
#' ecaim <- enhance_caim(caim, m)
#' plot(ecaim)
#' }
enhance_caim <- function(caim,
                         m,
                         w_red = 0.5,
                         sky_blue = sRGB(matrix(c(0, 0, 1), ncol = 3)),
                         thr = NULL,
                         fuzziness = NULL) {
  .check_if_r_was_normalized(caim, "caim")
  if (!compareRaster(caim, m, stopiffalse = FALSE)) {
    stop("\"x\" should match pixel by pixel whit \"m\".")
  }

  mem_sky_blue <- membership_to_color(caim, sky_blue)
  re_br <- caim$Red * w_red + caim$Blue * (1 - w_red)
  mem_sky_blue
  mem_thr <- suppressWarnings(local_fuzzy_thresholding(re_br,
                                                m,
                                                mem_sky_blue$membership_to_grey,
                                                thr,
                                                fuzziness))
  mem_thr[is.na(mem_thr)] <- 0
  mem <- mem_sky_blue$membership_to_target_color + mem_thr
  mem <- normalize(mem)
  names(mem) <- "Enhanced canopy image"
  mem
}
