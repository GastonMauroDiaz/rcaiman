#' Enhance canopy image
#'
#' This function was first proposed in \insertCite{Diaz2015;textual}{rcaiman}.
#' It uses the color perceptual attributes (hue, lightness, and chroma) to
#' enhance the contrast between the sky and plants through fuzzy classification.
#' It performs the next classification rules, here expressed in natural
#' language: clear sky is blue and clouds decrease its chroma; if clouds are
#' highly dense, then the sky is achromatic, and, in such cases, it can be light
#' or dark; everything that does not match this description is not sky. These
#' linguistic rules were translated to math language by means of fuzzy logic.
#'
#' This is a pixel-wise methodology that evaluates the possibility for a pixel
#' to be member of the class \emph{Gap}. High score could mean either high
#' membership to \code{sky_blue} or, in the case of achromatic pixels, a high
#' membership to values above \code{thr}. The algorithm internally uses
#' \code{\link{membership_to_color}} and \code{\link{local_fuzzy_thresholding}}.
#' The argument \code{sky_blue} is the \code{target_color} of the former
#' function, which output is the argument \code{mem} of the latter function.
#'
#' The argument \code{sky_blue} can be obtained from a photograph that clearly
#' shows the sky. Then, it can be used to process all the others taken with the
#' same equipment, configuration, and protocol.
#'
#' The \code{gamma} argument, along with \code{\link{gbc}}, is used to
#' back-correct the values passed to \code{\link{local_fuzzy_thresholding}}.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package.
#'
#' @inheritParams expand_noncircular
#' @inheritParams local_fuzzy_thresholding
#' @param m \linkS4class{SpatRaster}. A mask. For hemispherical photographs,
#'   check \code{\link{mask_hs}}. Default (\code{NULL}) is the equivalent to
#'   enter \code{!is.na(caim$Red)}. See the Details section in
#'   \code{\link{local_fuzzy_thresholding}} to understand how this argument can
#'   modify the output.
#' @param w_red Numeric vector of length one. Weight of the red channel. A
#'   single layer image is calculated as a weighted average of the blue and red
#'   channels. This layer is used as lightness information. The weight of the
#'   blue is the complement of \code{w_red}.
#' @param gamma Numeric vector of length one. This is for applying a gamma back
#'   correction to the lightness information (see Details and argument
#'   \code{w_red}).
#' @param sky_blue \linkS4class{color}. Is the \code{target_color} argument to
#'   be passed to \code{\link{membership_to_color}}. Default (\code{NULL}) is
#'   the equivalent to enter \code{sRGB(0.1, 0.4, 0.8)}--the HEX color code is
#'   #1A66CC, it can be entered into a search engine (such as Mozilla Firefox)
#'   to see a color swatch.
#'
#' @export
#' @references \insertAllCited{}
#' @family Pre-processing Functions
#' @return An object of class \linkS4class{SpatRaster}--with same pixel
#'   dimensions than \code{caim}--that should show more contrast between the sky
#'   and plant pixels than any of the individual bands from \code{caim}; if not,
#'   different parameters should be tested.
#'
#' @examples
#' \dontrun{
#' #circular hemispherical photo
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(1490, lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' blue <- caim$Blue %>% gbc()
#' plot(caim)
#'
#' sky_blue_sample <- read_caim(path, c(1092,1243), 66, 48)
#' plot(sky_blue_sample)
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>% normalize(.,0,255) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' hex(sky_blue)
#' # Use hex to obtain the HEX color code and enter it into a search
#' # engine (such as Mozilla Firefox) to see a color swatch. If the color is
#' # too pale (unsaturated), such as the one from the example (#6D90D0), it
#' # would be better to use the default. Alternatively, the values can be
#' # stretched, which often produces a more intense color. That is demonstrated
#' # below.
#' sky_blue_sample <- read_caim(path, c(1092,1243), 66, 48)
#' plot(sky_blue_sample)
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>% normalize() %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' hex(sky_blue) #005AFF
#'
#' caim <- normalize(caim)
#' ecaim <- enhance_caim(caim, m)
#' plot(ecaim)
#' plot(blue)
#'
#' m2 <- !mask_sunlit_canopy(caim, m) & m
#' hist(ecaim[m2])
#' hist(blue[m])
#'
#' plot(apply_thr(ecaim, thr_isodata(ecaim[m2])))
#' plot(apply_thr(blue, thr_isodata(blue[m])))
#'
#' #hemispherical photo from a smartphone
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' z <- zenith_image(2132/2, lens("Olloclip"))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)/2
#' caim <- expand_noncircular(caim, z, zenith_colrow) %>% normalize()
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- NA
#' blue <- caim$Blue %>% gbc()
#'
#' ecaim <- enhance_caim(caim, m)
#' plot(ecaim)
#' plot(blue)
#'
#' m2 <- !mask_sunlit_canopy(caim, m) & m
#' hist(ecaim[m2])
#' hist(blue[m])
#'
#' plot(apply_thr(ecaim, thr_isodata(ecaim[m2])))
#' plot(apply_thr(blue, thr_isodata(blue[m])))
#'
#' #restricted view canopy photo
#' path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#' caim <- read_caim(path)
#' plot(caim)
#' blue <- gbc(caim$Blue)
#' plot(blue)
#'
#' caim <- normalize(caim)
#' ecaim <- enhance_caim(caim)
#' plot(ecaim)
#'
#' m <- !mask_sunlit_canopy(caim)
#' hist(ecaim[])
#' hist(ecaim[m])
#' hist(blue)
#' plot(apply_thr(ecaim, thr_isodata(ecaim[m])))
#' plot(apply_thr(blue, thr_isodata(blue[])))
#'
#' }
enhance_caim <- function(caim,
                         m = NULL,
                         sky_blue = NULL,
                         w_red = 0,
                         thr = NULL,
                         fuzziness = NULL,
                         gamma = 2.2) {
  stopifnot(class(caim) == "SpatRaster")
  .was_normalized(caim, "caim")
  stopifnot(all(names(caim) == c("Red", "Green", "Blue")))
  if (is.null(m)) {
    m <- !is.na(caim$Red)
  } else {
    .is_single_layer_raster(m, "m")
    .is_logic_and_NA_free(m, "m")
    terra::compareGeom(caim, m)
  }


  if (is.null(sky_blue)) sky_blue <- colorspace::sRGB(0.1, 0.4, 0.8)

  mem_sky_blue <- membership_to_color(caim, sky_blue)
  lightness <- caim$Red * w_red + caim$Blue * (1 - w_red)
  if (!is.null(gamma)) lightness <- gbc(lightness * 255, gamma = gamma)
  mem_sky_blue
  mem_thr <- suppressWarnings(local_fuzzy_thresholding(lightness,
                                                m,
                                                mem_sky_blue$membership_to_grey,
                                                thr = thr,
                                                fuzziness = fuzziness))
  mem <- (mem_sky_blue$membership_to_target_color + mem_thr) / 2
  names(mem) <- "Enhanced canopy image"
  mem
}
