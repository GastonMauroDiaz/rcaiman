#' Out-of-the-box object-based image analysis of canopy photographs
#'
#' Out-of-the-box version of methods first presented in
#' \insertCite{Diaz2015;textual}{rcaiman}.
#'
#' This function is a hard-coded version of a pipeline that combines these main
#' functions \code{\link{mask_sunlit_canopy}}, \code{\link{enhance_caim}},
#' \code{\link{polar_qtree}}/\code{\link{qtree}}, and \code{\link{obia}}. The
#' code can be easily inspected by calling \code{ootb_obia} --no parenthesis.
#' Advanced users can use that code as a template.
#'
#' Pixels from the synthetic layer returned by \code{\link{obia}} that lay
#' between \code{0} and \code{1} are assigned to the class \emph{plant} only if
#' they are:
#'
#' \itemize{
#'
#' \item \code{0} after \code{\link{defuzzify}} with a sky grid segmentation of
#' \code{10} degrees.
#'
#' \item \code{0} after \code{\link{apply_thr}} with a threshold computed with
#' \code{\link{thr_isodata}}.
#'
#' \item Not exclusively surrounded by sky pixels.
#'
#' }
#'
#' Default values of \code{z} and \code{a} allows the processing of restricted
#' view photographs.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package.
#'
#' @inheritParams enhance_caim
#' @inheritParams ootb_mblt
#' @param m \linkS4class{SpatRaster}. Default (\code{NULL}) is the equivalent to
#'   enter \code{!is.na(z)} for hemispherical photography, or enter
#'   \code{!is.na(caim$Red)} for restricted view photography.
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#'
#'
#' @family Binarization Functions
#' @export
#'
#' @examples
#' \dontrun{
#' #circular hemispherical photo
#' path <- system.file("external/b4_2_5724.jpg", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2) %>%
#'   normalize()
#' z <- zenith_image(1490, lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#'
#' bin <- ootb_obia(caim, z, a)
#' plot(bin)
#'
#' ## to compare
#' blue <- gbc(caim$Blue*255)
#' plot(apply_thr(blue, thr_isodata(blue[!is.na(z)])))
#' plot(blue, col = seq(0,1,1/255) %>% grey())
#'
#' #hemispherical photo from a smartphone
#' path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#' caim <- read_caim(path) %>% normalize()
#' z <- zenith_image(2132/2, lens("Olloclip"))
#' a <- azimuth_image(z)
#' zenith_colrow <- c(1063, 771)/2
#' caim <- expand_noncircular(caim, z, zenith_colrow) %>% normalize()
#' m <- !is.na(caim$Red) & !is.na(z)
#' caim[!m] <- 0
#'
#' bin <- ootb_obia(caim, z, a)
#' plot(bin)
#'
#' ## to compare
#' blue <- gbc(caim$Blue*255)
#' plot(apply_thr(blue, thr_isodata(blue[m])))
#' plot(blue, col = seq(0,1,1/255) %>% grey())
#'
#' #restricted view canopy photo
#' path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#' caim <- read_caim(path) %>% normalize()
#'
#' bin <- ootb_obia(caim)
#' plot(bin)
#'
#' ## to compare
#' blue <- gbc(caim$Blue*255)
#' plot(apply_thr(blue, thr_isodata(blue[])))
#' plot(blue, col = seq(0,1,1/255) %>% grey())
#' }
ootb_obia <- function(caim, z = NULL, a = NULL, m = NULL, sky_blue = NULL) {

  if (is.null(m)) {
    if (is.null(z)) {
      m <- !is.na(caim$Red)
    } else {
      m <- !is.na(z)
    }
  }

  m2 <- !mask_sunlit_canopy(caim, m)
  ecaim <- enhance_caim(caim, m, sky_blue = sky_blue,
                        w_red = 0, gamma = 2.2, thr = NULL,
                        fuzziness = NULL)
  bin <- apply_thr(ecaim, thr_isodata(ecaim[m2]))

  if (is.null(z)) {
    seg <- qtree(caim, scale_parameter = 0.2)
    size <- round(max(c(ncol(caim), nrow(caim)))/22)
    g <- chessboard(caim, size)
  } else {
    if (is.null(a)) a <- azimuth_image(z)
    seg <- polar_qtree(caim, z, a, scale_parameter = 0.2)
    g <- sky_grid_segmentation(z, a, 10)
  }

  r <- gbc(caim$Blue*255, gamma = 2.2)
  synth <- obia(r, z, a, bin, seg)
  foliage <- !is.na(synth)
  synth <- terra::cover(synth, bin)

  bin_obia <- defuzzify(synth, g) | apply_thr(synth, thr_isodata(synth[foliage]))
  ma <- matrix(c(1,1,1,1,-8,1,1,1,1), ncol = 3, nrow = 3)
  bin_obia[terra::focal(bin_obia, ma) == 8] <- 1
  bin_obia[!bin] <- 0
  as.logical(bin_obia)
}
