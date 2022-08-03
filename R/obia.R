#' Object-based image analysis of hemispherical photographs
#'
#' This method was first presented in \insertCite{Diaz2015;textual}{rcaiman}.
#' This version is simpler since it relies on a better working binarized image.
#' The version from 2015 uses an automatic selection of samples followed by a
#' \emph{knn} classification of segments containing foliage. This version uses
#' de gap fraction extracted from \code{bin} to classify \emph{foliage} by
#' defining upper and lower limits through the arguments \code{gf_mx} and
#' \code{gf_mn}.
#'
#' This method produces a synthetic layer by computing the ratio of \code{r} to
#' the maximum value of \code{r} at the segment level. This process is carried
#' out only on the pixels covered by the classes \emph{foliage} and \emph{sky};
#' the latter is defined by bin equal to one. To avoid spurious values, the
#' quantile \code{0.9} is computed instead of the maximum. Pixels from the
#' synthetic layer between \code{0} and \code{1} are binarized by following two
#' criteria, i.e., to assign a pixel to the class \emph{plant} it has to be
#' \code{0} under both criteria. Those criteria are (1) \code{defuzzify} with a
#' sky grid segmentation of \code{10} degrees and (2) \code{apply_thr} with a
#' threshold of \code{0.5}. In addition, to be \emph{plant}, the pixel cannot be
#' exclusively surrounded by sky pixels.
#'
#' @inheritParams ootb_mblt
#' @param segmentation \linkS4class{SpatRaster} built with
#'   \code{\link{polar_qtree}}.
#' @param gf_mn,gf_mx Numeric vector of length one. The minimum/maximum gap
#'   fraction that a segment should comply with to be considered as one
#'   containing foliage.
#'
#' @return \linkS4class{SpatRaster}. An improved version of \code{bin}.
#' @export
#'
#' @family Binarization Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' sky_blue_sample <- crop(caim, ext(610,643,760,806))
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' ecaim <- enhance_caim(caim, !is.na(z), sky_blue, gamma = 2.2)
#' bin <- apply_thr(ecaim, 0.5)
#' seg <- polar_qtree(caim, z, a)
#' bin_obia <- obia(gbc(caim$Blue), z, a, bin, seg)
#' plot(bin - bin_obia)
#' }
obia <- function(r, z, a, bin, segmentation, gf_mn = 0.2, gf_mx = 0.95) {
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  gf <- extract_feature(bin, segmentation)
  gf1 <- extract_feature(bin, sky_grid_segmentation(z, a, 1))
  foliage <- (gf > gf_mn & gf < gf_mx)
  foliage[gf1 == 0] <- 0
  segmentation <- segmentation * foliage
  segmentation[segmentation == 0] <- NA
  synth <- r/extract_feature(r, segmentation, function(x) quantile(x, 0.9))
  synth[gf == 1] <- 1
  synth[synth > 1 | synth < 0] <- NA
  synth[is.na(synth)] <- 1
  synth[!bin] <- 0

  g <- sky_grid_segmentation(z, a, 10)
  bin_obia <- defuzzify(synth, g) | apply_thr(synth, 0.5)
  ma <- matrix(c(1,1,1,1,-8,1,1,1,1), ncol = 3, nrow = 3)
  bin_obia[terra::focal(bin_obia, ma) == 8] <- 1
  bin_obia[!bin] <- 0
  as.logical(bin_obia)
}
