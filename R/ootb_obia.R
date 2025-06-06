#' Out-of-the-box object-based image analysis of canopy photographs
#'
#' Out-of-the-box version of methods first presented in
#' \insertCite{Diaz2015;textual}{rcaiman}.
#'
#' This function is a hard-coded version of a pipeline that combines these main
#' functions [mask_sunlit_canopy()], [enhance_caim()],
#' [polar_qtree()]/[qtree()], and [obia()]. The code can be easily inspected by
#' calling `ootb_obia` --no parenthesis. Advanced users can use that code as a
#' template.
#'
#' Pixels from the synthetic layer returned by [obia()] that lay between `0` and
#' `1` are assigned to the class *plant* only if they comply with the following
#' conditions:
#'
#' * Their values are equal to `0` after [defuzzify()] with a
#' sky grid segmentation of `10` degrees.
#' * Their values are equal to `0` after [apply_thr()] with a
#' threshold computed with [thr_isodata()].
#' * They are not exclusively surrounded by sky pixels.
#'
#' Use the default values of `z` and `a` to process restricted view photographs.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} or
#' \insertCite{Diaz2023;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`)).
#'
#' @inheritParams enhance_caim
#' @inheritParams obia
#' @param m [SpatRaster-class]. Default (`NULL`) is the equivalent to enter
#'   `!is.na(z)` for hemispherical photography, or enter `!is.na(caim$Red)` for
#'   restricted view photography.
#'
#' @return An object of class [SpatRaster-class] with values `0` and `1`.
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#'
#' # =====================================
#' # Hemispherical Photo from a Smartphone
#' # =====================================
#'
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path) %>% normalize_minmax()
#' z <- zenith_image(2132/2, c(0.7836, 0.1512, -0.1558))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)/2
#' caim <- expand_noncircular(caim, z, zenith_colrow) %>% normalize_minmax()
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#'
#' bin <- ootb_obia(caim, z, a)
#' plot(bin)
#'
#' # ============================
#' # Restricted View Canopy Photo
#' # ============================
#'
#' path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#' caim <- read_caim(path) %>% normalize_minmax()
#'
#' bin <- ootb_obia(caim)
#' plot(bin)
#' }
ootb_obia <- function(caim, z = NULL, a = NULL, m = NULL,
                      sky_blue = NULL, w_red = 0, gamma = 2.2) {

  if (is.null(m)) {
    if (is.null(z)) {
      m <- !is.na(caim$Red)
    } else {
      m <- !is.na(z)
    }
  }

  ecaim <- enhance_caim(caim, m, sky_blue = sky_blue,
                        w_red = w_red, gamma = gamma, thr = NULL,
                        fuzziness = NULL)
  if (is.null(gamma)) {
    bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
  } else {
    m2 <- !mask_sunlit_canopy(caim, m) & m
    bin <- apply_thr(ecaim, thr_isodata(ecaim[m2]))
  }


  if (is.null(z)) {
    seg <- qtree(caim, scale_parameter = 0.2)
    size <- round(max(c(ncol(caim), nrow(caim)))/22)
    g <- chessboard(caim, size)
  } else {
    if (is.null(a)) a <- azimuth_image(z)
    seg <- polar_qtree(caim, z, a, scale_parameter = 0.2)
    g <- sky_grid_segmentation(z, a, 10)
  }

  if (is.null(gamma)) {
    bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
  } else {
    bin <- apply_thr(ecaim, thr_isodata(ecaim[m2]))
  }

  if (is.null(gamma)) {
    r <- caim$Blue
  } else {
    r <- gbc(caim$Blue*255, gamma = gamma)
  }

  synth <- obia(r, z, a, bin, seg)
  foliage <- !is.na(synth)
  synth <- terra::cover(synth, bin)

  bin_obia <- defuzzify(synth, g) | apply_thr(synth, thr_isodata(synth[foliage]))
  rem_isolated_black_pixels(bin_obia)
}
