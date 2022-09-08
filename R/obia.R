#' Object-based image analysis of canopy photographs
#'
#' Object-based image analysis targeting the canopy  silhouette.
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
#' out only on the pixels covered by the classes \emph{foliage} and \emph{sky}--
#' the latter is defined by \code{bin} equal to one. To avoid spurious values,
#' the quantile \code{0.9} is computed instead of the maximum. Pixels not
#' belonging to the class \emph{foliage} return as \code{NA}.
#'
#' Default values of \code{z} and \code{a} allows the processing of
#' restricted view photographs.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package.
#'
#' @inheritParams ootb_mblt
#' @param bin  \linkS4class{SpatRaster}. This should be a working binarization
#'   of \code{r} without gross errors.
#' @param segmentation \linkS4class{SpatRaster} built with
#'   \code{\link{polar_qtree}} or \code{\link{qtree}}.
#' @param gf_mn,gf_mx Numeric vector of length one. The minimum/maximum gap
#'   fraction that a segment should comply with to be considered as one
#'   containing foliage.
#'
#' @return \linkS4class{SpatRaster}.
#' @export
#'
#' @family Binarization Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue %>% gbc()
#' caim <- normalize(caim)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' m2 <- !mask_sunlit_canopy(caim, m)
#' ecaim <- enhance_caim(caim, m)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m2]))
#'
#' seg <- polar_qtree(caim, z, a)
#' synth <- obia(r, z, a, bin, seg)
#' plot(synth)
#' foliage <- !is.na(synth)
#' hist(synth[foliage])
#' synth <- terra::cover(synth, bin)
#' plot(synth)
#' bin_obia <- apply_thr(synth, thr_isodata(synth[foliage]))
#' plot(bin - bin_obia)
#' plot(bin_obia)
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
