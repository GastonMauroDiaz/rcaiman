#' Object-based image analysis of canopy photographs
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
#' the latter is defined by bin equal to one. To avoid spurious values, the
#' quantile \code{0.9} is computed instead of the maximum. Pixels not belonging
#' to the class \emph{foliage} return as \code{NA}.
#'
#' @inheritParams ootb_mblt
#' @param bin  \linkS4class{SpatRaster}. This should be a working binarization
#'   of \code{r} without gross errors.
#' @param segmentation \linkS4class{SpatRaster} built with
#'   \code{\link{polar_qtree}}.
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
#' seg <- polar_qtree(caim, z, a)
#' ecaim <- enhance_caim(caim, !is.na(z))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[]))
#' synth <- obia(r, z, a, bin, seg)
#' plot(synth)
#' synth <- terra::cover(synth, bin)
#' plot(synth)
#' synth[!bin] <- 0
#' plot(synth)
#' bin_obia <- apply_thr(synth, thr_isodata(synth[]))
#' plot(bin - bin_obia)
#' plot(bin_obia)
#' # This is the closest to the workflow presented in Diaz and lencinas
#' # (2015) that this package allows.
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
  # synth[gf == 1] <- 1
  synth[synth > 1 | synth < 0] <- NA
  synth[is.na(synth)] <- 1
  # synth[!bin] <- 0
  synth[!foliage] <- NA
  synth
}
