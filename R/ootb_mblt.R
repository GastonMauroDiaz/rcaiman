#' Out-of-the-box model-based local thresholding
#'
#' Out-of-the-box version of the model-based local thresholding (MBLT)
#' algorithm.
#'
#' This function is a hard-coded version of a MBLT pipeline that starts
#' producing a working binarized image and ends with a refined binarized image.
#' The pipeline combines these main functions \code{\link{find_sky_pixels}}--if
#' \code{bin} is \code{NULL}--, \code{\link{fit_coneshaped_model}},
#' \code{\link{find_sky_pixels_nonnull}}, and \code{\link{fit_trend_surface}}.
#' The code can be easily inspected by calling \code{ootb_mblt}--no
#' parenthesis. Advanced users can use that code as a template.
#'
#' The MBLT algorithm was first presented in
#' \insertCite{Diaz2018;textual}{rcaiman}. The version presented here differs
#' from that in the following main aspects:
#'
#' \itemize{
#'
#' \item \code{intercept} is set to \code{0}, \code{slope} to \code{1}, and
#' \code{w} to \code{0.5}
#'
#' \item This version implements a regional thresholding approach as the first
#' step instead of a global one. Please refer to \code{\link{find_sky_pixels}}.
#'
#' \item It does not use asynchronous acquisition under the open sky. So, the
#' cone-shaped model (\code{\link{fit_coneshaped_model}}) run without a filling
#' source, but the result of it is used as filling source for trend surface
#' fitting (\code{\link{fit_trend_surface}}).
#'
#' \item  \code{\link{find_sky_pixels_nonnull}} is used to update the first
#' working binarized image, after \code{\link{fit_coneshaped_model}}.
#'
#' }
#'
#' This function searches for black objects against a light background. When
#' regular canopy hemispherical images are provided as input, the algorithm will
#' find dark canopy elements against a bright sky almost everywhere in the
#' picture and, therefore, the result will fit user expectations. However, if a
#' hemispherical photograph taken under the open sky is provided, this algorithm
#' would be still searching black objects against a light background, so the
#' darker portions of the sky will be taken as objects, i.e., canopy. As a
#' consequence, this will not fit users expectations since they are looking for
#' the classes \emph{Gap} and \emph{No-gap}, no matter if one of those are not
#' in the picture itself. This kind of error could happen with photographs of
#' open forests for the same working principle.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018;textual}{rcaiman} in addition to this package.
#'
#' @param r \linkS4class{SpatRaster}. A normalized greyscale image. Typically,
#'   the blue channel extracted from a canopy photograph. Please see
#'   \code{\link{read_caim}} and \code{\link{normalize}}.
#' @param z \linkS4class{SpatRaster} built with \code{\link{zenith_image}}.
#' @param a \linkS4class{SpatRaster} built with \code{\link{azimuth_image}}.
#' @param bin \linkS4class{SpatRaster}. This should be a preliminary
#'   binarization of \code{r} useful for masking pixels that are very likely
#'   to be pure sky pixels.
#'
#' @export
#' @family Binarization Functions
#' @seealso \code{\link{thr_image}}
#'
#' @return Object from class list containing the binarized image (named
#'   \emph{bin}) and the reconstructed skies (named \emph{sky_cs} and
#'   \emph{sky_s}).
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' r[is.na(z)] <- 0 #because FOV > 180
#' bin <- ootb_mblt(r, z, a)
#' plot(bin$bin)
#' }
ootb_mblt <- function(r, z, a, bin = NULL) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  if (is.null(bin)) {
    bin <- find_sky_pixels(r, z, a)
  }
  g <- sky_grid_segmentation(z, a, 10)
  sky_points <- extract_sky_points(r, bin, g)
  sky_points <- extract_rl(r, z, a, sky_points, NULL)
  model <- fit_coneshaped_model(sky_points$sky_points)
  if (is.null(model)) {
    sky_cs <- z
    terra::values(sky_cs) <- median(sky_points$sky_points$dn)
  } else {
    sky_cs <- model$fun(z, a)
    bin <- find_sky_pixels_nonnull(r, sky_cs, g)
  }

  sky_s <- fit_trend_surface(r, z, a, bin,
                             filling_source = sky_cs,
                             np = 6)$image
  thr <- suppressWarnings(thr_image(sky_s, 0, 0.5))
  bin <- apply_thr(r, thr)
  list(bin = bin, sky_cs = sky_cs, sky_s = sky_s)
}
