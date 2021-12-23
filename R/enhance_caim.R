#' Enhance canopy image
#'
#' This function is presented in \insertCite{Diaz2015;textual}{rcaiman}. It uses
#' the color perceptual attributes to enhance the contrast between the sky and
#' plants through fuzzy classification. Color has three different perceptual
#' attributes: hue, lightness, and chroma. The algorithm was developed using the
#' following premise: the color of the sky is different from the color of
#' plants. It performs the next classification rules, here expressed in natural
#' language: clear sky is blue and clouds decrease its chroma; if clouds are
#' highly dense, then the sky is achromatic, and, in such cases, it can be light
#' or dark; everything that does not match this description is not sky. These
#' linguistic rules were translated to math language by means of fuzzy logic.
#'
#' @inheritParams expand_noncircular
#' @inheritParams local_fuzzy_thresholding
#' @param w_red Numeric vector of length one. Weight of the red channel. A
#'   single layer image is calculated as a weighted average of the blue and red
#'   channels. This layer is used as lightness information. The weight of the
#'   blue is the complement of \code{w_red}.
#' @param sky_blue \linkS4class{color}. Is the \code{target_color} argument to
#'   be passed to \code{\link{membership_to_color}}. By default is pure Blue.
#'
#' @details This is a pixel-wise methodology that evaluates the possibility for
#'   a pixel to be member of the class Gap. High score could mean either high
#'   membership to \code{sky_blue} or, in the case of achromatic pixels, a high
#'   membership to \code{thr}. The algorithm internally uses
#'   \code{\link{membership_to_color}} and
#'   \code{\link{local_fuzzy_thresholding}}. The argument \code{sky_blue} is the
#'   \code{target_color} of the former function, which output is the argument
#'   \code{mem} of the latter function.
#'
#' @export
#' @references \insertRef{Diaz2015}{rcaiman}
#' @family Fuzzy logic functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' m <- !is.na(z)
#' ecaim <- enhance_caim(caim, m)
#' plot(ecaim)
#' }
enhance_caim <- function(caim,
                         m,
                         w_red = 0.5,
                         sky_blue = NULL) {
  .check_if_r_was_normalized(caim, "caim")
  if (!compareRaster(caim, m, stopiffalse = FALSE)) {
    stop("\"x\" should match pixel by pixel whit \"m\".")
  }

  if (is.null(sky_blue)) {
      sky_blue <- colorspace::sRGB(matrix(c(0, 0, 1), ncol = 3))
    }
  mem_sky_blue <- membership_to_color(caim, sky_blue)
  re_br <- caim$Red * w_red + caim$Blue * (1 - w_red)
  mem_sky_blue
  mem_thr <- suppressWarnings(local_fuzzy_thresholding(re_br,
                                                m,
                                                mem_sky_blue$membership_to_grey,
                                                thr = NULL,
                                                fuzziness = NULL))
  mem_thr[is.na(mem_thr)] <- 0
  mem <- mem_sky_blue$membership_to_target_color + mem_thr
  mem <- normalize(mem)
  names(mem) <- "Enhanced canopy image"
  mem
}
