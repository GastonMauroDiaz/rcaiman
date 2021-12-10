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

#' Compute the membership value to a color
#'
#' Compute the degree of membership to a color with two Gaussian membership
#' functions and axis \emph{A} and \emph{B} from the \emph{CIE L*a*b*} color
#' space. The lightness information is not considered in the calculations.
#'
#' @param r \linkS4class{RasterStackBrick}.
#' @param target_color \linkS4class{color}.
#' @param sigma Numeric vector of length one. Use \code{NULL} (default) to order
#'   to estimate it automatically as the euclidean distance between
#'   \code{target_color} and grey in the \emph{CIE L*a*b*} color space.
#'
#' @export
#'
#' @examples
membership_to_color <- function(r, target_color, sigma = NULL) {
  .check_if_r_was_normalized(r)
  stopifnot(names(r) == c("Red", "Green", "Blue"))
  color <- colorspace::sRGB(values(r))
  error_msn <- paste("\"target_color\" must be a subclass of the",
                     "virtual class named \"color\" from",
                     "\"colorspace\" package")
  .is_class_from_colorspace <- function(x) {
    if (!isS4(x)) stop(error_msn)
    if (attr(class(x), "package") != "colorspace") stop(error_msn)
  }
  .is_class_from_colorspace(color)
  .is_class_from_colorspace(target_color)
  if (class(x) != "LAB") x <- as(x, "LAB")
  p <- .get_gaussian_2d_parameters(target_color, sigma)
  max_z <- .gaussian2d(p[1], p[2], p[1], p[2], p[3])
  x <- colorspace::coords(x)
  m <- .gaussian2d(x[, 2], x[, 3], p[1], p[2], p[3]) / max_z
  r$Blue <- unname(m)
  r <- r$Blue
  names(r) <- "Membership to the target color"
  r
}

#' Compute a weighted membership value to lightness
#'
#' This function uses a threshold value as the location parameter of a logistic
#' membership function whose scale parameter depends on a variable. This
#' dependence can be explained as follows: if the variable is equal to \code{1},
#' then the membership function is as a threshold function because the scale
#' parameter is \code{0}; lowering the variable increases the scale parameter,
#' thus blurring the threshold because it decreases the steepness of the curve.
#'
#' @param lightness Numeric vector. The lightness value.
#' @param m Numeric vector lying between \code{0} and \code{1}, same length as
#'   \code{lightness}. It is the scale parameter of the logistic membership
#'   function.
#' @param thr Numeric vector of length one. Location parameter of the logistic
#'   membership function.
#' @param fuzziness Numeric vector of length one.
#'
#' @noRd
fuzzy_lightness <- function (lightness, m, thr, fuzziness) {
  stats::plogis(lightness, thr, fuzziness * (1 - m))
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
#' @param r
#' @param m
#' @param w_red
#' @param w_blue
#' @param sky_blue
#' @param thr
#' @param fuzziness
#'
#' @details This is a pixel-wise methodology that evaluates if pixels are sky
#'   blue colored. High score means high membership to \code{sky_blue}. When a
#'   pixel is achromatic, then it uses pixel brightness. The algorithm
#'   internally uses \code{\link{membership_to_color}} and
#'   \code{\link{fuzzy_lightness}}. The argument \code{sky_blue} is the
#'   \code{target_color} of the former function, which output is the argument
#'   \code{m} of the latter function. To evaluate the brightness of an
#'   achromatic pixel, the algorithm uses \strong{Relative Brightness} (see
#'   references).
#'
#'   Argument \code{mask} can be used to affect the estimation of two arguments
#'   of \code{\link{fuzzyLightness}}. Affected arguments are \code{thr} and
#'   \code{fuzziness}. The function \code{\link{autoThr}} is used to estimate
#'   \code{thr}. To compute \code{fuzziness}, the algorithm takes the maximum
#'   and the minimum values of the Relative Brightness and calculate its mean.
#'
#' @return
#' @export
#'
#' @examples
enhance_caim <- function(r,
                         m,
                         w_red = 0.25,
                         w_blue = 0.75,
                         sky_blue = sRGB(matrix(c(0.529, 0.808, 0.921), ncol = 3)),
                         thr = NULL,
                         fuzziness = NULL) {
  .check_if_r_was_normalized(r)
  if (!compareRaster(r, m, stopiffalse = FALSE)) {
    stop("\"x\" should match pixel by pixel whit \"m\".")
  }

  m_sky_blue <- membership_to_color(colorspace::sRGB(as.matrix(r)), sky_blue)

  .relative_brightness <- function(x, w_red, w_blue) {

  }

  re_br <- r$Red * w_red + r$Blue * w_blue

  if (is.null(thr)) {
    if (!requireNamespace("autothresholdr", quietly = TRUE)) {
      stop(paste(
        "Package \"autothresholdr\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
      )
    }
    dns <- re_br[m]
    thr <- autothresholdr::auto_thresh(round(dns * 255), "IsoData")[1] / 255
  }
  if (is.null(fuzziness)) {
    fuzziness <- (max(re_br[m]) - min(re_br[m])) / 2
  }

  m_light <- fuzzy_lightness(values(re_br), m_sky_blue, thr, fuzziness)
  membership <- m_sky_blue * m_light
  r$Blue[] <- membership
  r <- r$Blue
  names(r) <- "Enhanced canopy image"
  r
}



