#' Enhance canopy image
#'
#' This function was first proposed in \insertCite{Diaz2015;textual}{rcaiman}.
#' It uses the color perceptual attributes (hue, lightness, and chroma) to
#' enhance the contrast between the sky and plants through fuzzy classification.
#' It applies the next classification rules (here expressed in natural
#' language): clear sky is blue and clouds decrease its chroma; if clouds are
#' highly dense, then the sky is achromatic, and, in such cases, it can be light
#' or dark; everything that does not match this description is not sky. These
#' linguistic rules were translated to math language by means of fuzzy logic.
#' This translation was thoughtfully explained in the aforementioned article.
#'
#' This is a pixel-wise methodology that evaluates the possibility for a pixel
#' to be member of the class *Gap*. High score could mean either high membership
#' to `sky_blue` or, in the case of achromatic pixels, a high membership to
#' values above `thr`. The algorithm internally uses [membership_to_color()] and
#' [local_fuzzy_thresholding()]. The argument `sky_blue` is the `target_color`
#' of the former function and its output is the argument `mem` of the latter
#' function.
#'
#' The argument `sky_blue` can be obtained from a photograph that clearly shows
#' the sky. Then, it can be used to process all the others photograph taken with
#' the same equipment, configuration, and protocol.
#'
#' Via the `gamma` argument, [gbc()] can be internally used to back-correct the
#' values passed to [local_fuzzy_thresholding()].
#'
#' @note If you use this function in your research, please cite
#'   \insertCite{Diaz2015;textual}{rcaiman} in addition to this package
#'   (`citation("rcaiman"`).
#'
#' @inheritParams expand_noncircular
#' @inheritParams local_fuzzy_thresholding
#' @inheritParams membership_to_color
#' @param m [SpatRaster-class]. A mask. For hemispherical photographs, check
#'   [select_sky_vault_region()]. Default (`NULL`) is the equivalent to enter
#'   `!is.na(caim$Red)`.
#' @param w_red Numeric vector of length one. Weight of the red channel. A
#'   single layer image is calculated as a weighted average of the blue and red
#'   channels. This layer is used as lightness information. The weight of the
#'   blue is the complement of `w_red`.
#' @param gamma Numeric vector of length one. This is for applying a gamma back
#'   correction to the lightness information (see Details and argument `w_red`).
#' @param sky_blue [color-class]. Is the `target_color` argument to be passed to
#'   [membership_to_color()]. Default (`NULL`) is the equivalent to enter
#'   `sRGB(0.1, 0.4, 0.8)`--the HEX color code is #1A66CC, it can be entered
#'   into a search engine (such as Mozilla Firefox) to see a color swatch.
#'
#' @note The default value of argument `m` is the equivalent to enter
#'   `!is.na(caim$Red)`. See the Details section in [local_fuzzy_thresholding()]
#'   to understand how this argument can modify the output.
#'
#' @export
#' @references \insertAllCited{}
#'
#' @return An object of class [SpatRaster-class] (with same pixel dimensions
#'   than `caim`) that should show more contrast between the sky and plant
#'   pixels than any of the individual bands from `caim`; if not, different
#'   parameters should be tested.
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r <- normalize_minmax(caim$Blue)
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' mx <- optim_max(caim, bin)
#' mn <- min(caim[m])
#'
#' sky_blue_sample <- crop_caim(caim, c(327, 239), 41, 89)
#' plotRGB(normalize_minmax(sky_blue_sample, mn, mx, TRUE)*255)
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>%
#'   normalize_minmax(., mn, mx) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' hex(sky_blue)
#' # Use hex() to obtain the HEX color code. To see a color swatch, enter the
#' # HEX code into a search engine (such as Mozilla Firefox).
#' # NOTE: see extract_dn() for a better method to obtain sky_blue
#' sky_blue <- polarLAB(50, 17, 293)
#'
#' caim <- normalize_minmax(caim, mx = mx, force_range = TRUE)
#' ecaim <- enhance_caim(caim, m, sky_blue = sky_blue)
#' plot(ecaim)
#' plot(caim$Blue)
#'
#' ## to compare
#' plot(apply_thr(ecaim, thr_isodata(ecaim[m])))
#' plot(apply_thr(caim$Blue, thr_isodata(caim$Blue[m])))
#' }
enhance_caim <- function(caim,
                         m = NULL,
                         sky_blue = NULL,
                         sigma = NULL,
                         w_red = 0,
                         thr = NULL,
                         fuzziness = NULL,
                         gamma = NULL) {
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

  if (is.null(sky_blue)) {
    sky_blue <- colorspace::sRGB(0.1, 0.4, 0.8)
  }

  mem_sky_blue <- membership_to_color(caim, sky_blue, sigma)
  lightness <- caim$Red * w_red + caim$Blue * (1 - w_red)
  if (!is.null(gamma)) {
    stopifnot(is.numeric(gamma))
    lightness <- gbc(lightness * 255, gamma = gamma)
  }
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
