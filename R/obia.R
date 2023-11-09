#' Do object-based image analysis of canopy photographs
#'
#' Object-based image analysis targeting the canopy  silhouette.
#'
#' This method was first presented in \insertCite{Diaz2015;textual}{rcaiman}.
#' This version is simpler since it relies on a better working binarized image.
#' The version from 2015 uses an automatic selection of samples followed by a
#' *knn* classification of segments containing foliage. This version uses
#' de gap fraction extracted from `bin` to classify *foliage* by defining upper
#' and lower limits through the arguments `gf_mx` and `gf_mn`.
#'
#' This method produces a synthetic layer by computing the ratio of `r` to the
#' maximum value of `r` at the segment level. This process is carried out only
#' on the pixels covered by the classes *foliage* and *sky*. The latter is
#' defined by `bin` equal to one. To avoid spurious values, the quantile `0.9`
#' is computed instead of the maximum. Pixels not belonging to the class
#' *foliage* return as `NA`.
#'
#' Default values of `z` and `a` allows the processing of restricted view
#' photographs.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`).
#'
#' @inheritParams ootb_mblt
#' @param bin  [SpatRaster-class]. This should be a working binarization of `r`
#'   without gross errors.
#' @param segmentation [SpatRaster-class] built with [polar_qtree()] or
#'   [qtree()].
#' @param gf_mn,gf_mx Numeric vector of length one. The minimum/maximum gap
#'   fraction that a segment should comply with to be considered as one
#'   containing foliage.
#'
#' @return [SpatRaster-class].
#' @export
#'
#' @family Binarization Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' ecaim <- enhance_caim(caim, m)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' plot(bin)
#'
#' seg <- polar_qtree(caim, z, a)
#' synth <- obia(caim$Blue, z, a, bin, seg)
#' plot(synth)
#' foliage <- !is.na(synth)
#' hist(synth[foliage])
#' synth <- terra::cover(synth, bin)
#' plot(synth)
#' hist(synth[foliage])
#' }
obia <- function(r, z = NULL, a = NULL,
                 bin, segmentation, gf_mn = 0.2, gf_mx = 0.95) {
  .is_single_layer_raster(r, "r")
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")

  gf <- extract_feature(bin, segmentation)
  foliage <- (gf > gf_mn & gf < gf_mx)
  segmentation <- segmentation * foliage
  segmentation[segmentation == 0] <- NA
  synth <- r/extract_feature(r, segmentation,
                             function(x) quantile(x, 0.9, na.rm = TRUE))
  synth[synth > 1 | synth < 0] <- NA
  synth[is.na(synth)] <- 1
  synth[!foliage] <- NA
  synth
}
