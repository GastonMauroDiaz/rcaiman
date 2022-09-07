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
#'
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
#' #hemispherical photo from a smartphone
#' caim <- read_caim()
#' r <- caim$Blue %>% gbc()
#' caim <- normalize(caim)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' seg <- polar_qtree(caim, z, a)
#' ecaim <- enhance_caim(caim, !is.na(z))
#' m <- mask_sunlit_canopy(caim, !is.na(z))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[!m]))
#' synth <- obia(r, z, a, bin, seg)
#' plot(synth)
#' foliage <- !is.na(synth)
#' synth <- terra::cover(synth, bin)
#' plot(synth)
#' synth[!bin] <- 0
#' plot(synth)
#' bin_obia <- apply_thr(synth, thr_isodata(synth[foliage]))
#' plot(bin - bin_obia)
#' plot(bin_obia)
#' ## This is the closest to the workflow presented in Diaz and Lencinas
#' ## (2015) that this package allows.
#'
#' #restricted view canopy photo
#' path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#'
#' sky_blue_sample <- read_caim(path, c(535,276), 11, 12) %>% normalize(., 0, 255)
#' plotRGB(sky_blue_sample*255)
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#'
#' caim <- read_caim(path)
#' plot(caim)
#' caim <- normalize(caim)
#' ecaim <- enhance_caim(caim, sky_blue = sky_blue)
#' m <- mask_sunlit_canopy(caim)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[!m]))
#'
#' seg <- qtree(caim)
#' synth <- obia(blue, NULL, NULL, bin, seg)
#' plot(synth)
#' foliage <- !is.na(synth)
#' synth <- cover(synth, bin)
#' plot(synth)
#' synth[!bin] <- 0
#' plot(synth)
#' bin_obia <- apply_thr(synth, thr_isodata(synth[foliage]))
#' plot(bin - bin_obia)
#' plot(bin_obia)
#' }
obia <- function(r, z = NULL, a = NULL, bin, segmentation, gf_mn = 0.2, gf_mx = 0.95) {
  .is_single_layer_raster(r, "r")
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")

  if (is.null(z)) {
    size <- round(max(c(ncol(r), nrow(r)))/250)
    if (size < 2) size <- 2
    g1 <- chessboard(r, size)
  } else {
    g1 <- sky_grid_segmentation(z, a, 1)
  }
  gf1 <- extract_feature(bin, g1)
  gf <- extract_feature(bin, segmentation)
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
