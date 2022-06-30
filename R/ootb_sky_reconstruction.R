#' Out-of-the-box sky reconstruction
#'
#' Build an above canopy image from a single below canopy image.
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions: \code{\link{ootb_mblt}}, \code{\link{fit_cie_sky_model}}, and
#' \code{\link{interpolate_sky_points}}. The code can be easily inspected by
#' calling \code{ootb_sky_reconstruction} --no parenthesis. Advanced users could
#' use that code as a template.
#'
#' This pipeline is based on these two studies:
#' \insertCite{Lang2010;textual}{rcaiman} and
#' \insertCite{Diaz2018;textual}{rcaiman}.
#'
#' Details about how the original method by
#' \insertCite{Diaz2018;textual}{rcaiman} was adapted can be read here
#' \code{\link{ootb_mblt}}.
#'
#' The main differences between the original method by
#' \insertCite{Lang2010;textual}{rcaiman} and the one implemented here are: (1)
#' it is fully automatic, (2) the residuals of the CIE sky model
#' ($residuals=model-data$) are interpolated instead of the sky digital numbers
#' (the data), and (3) the final sky reconstruction is obtained by subtracting
#' the interpolated residuals to the CIE sky model instead of by calculating a
#' weighted average parameterized by the user.
#'
#' The recommended input for this function is data pre-processed with the HSP
#' software package \insertCite{Lang2013}{rcaiman}. Please, refer to
#' code{link{write_sky_marks}} for additional details about HSP, and refer to
#' \code{\link{fit_cie_sky_model}} and \code{\link{interpolate_sky_points}} to
#' know why the HSP pre-processing is convenient.
#'
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_trend_surface
#'
#'
#' @export
#'
#' @references \insertRef{Lang2010}{rcaiman}
#'
#'   \insertRef{Lang2013}{rcaiman}
#'
#'   \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(caim$Blue)
#' sky <- ootb_sky_reconstruction(blue, z, a)
#' plot(sky)
#' ratio <- blue / sky
#' plot(ratio)
#' hist(ratio)
#' }
ootb_sky_reconstruction <- function(r, z, a, bin = NULL) {
  if (is.null(bin)) {
    bin <- ootb_mblt(r, z, a)$bin
  }
  bin[mask_hs(z, 80, 90)] <- 0
  bin <- as.logical(bin)
  g <- sky_grid_segmentation(z, a, 10)
  sky_points <- extract_sky_points(r, bin, g)
  zenith_dn <- extract_zenith_dn(r, z, a, sky_points)
  sun_coord <- extract_sun_coord(r, z, a, bin, g, max_angular_dist = 45)
  model <- fit_cie_sky_model(r, z, a,
                             zenith_dn$sky_points,
                             zenith_dn$zenith_dn,
                             sun_coord)
  sky_cie <- model$relative_luminance * model$zenith_dn

  residu <- sky_cie - r
  zenith_dn <- suppressWarnings(extract_zenith_dn(residu, z, a, sky_points))
  zenith_dn$sky_points$rl <- zenith_dn$sky_points$dn
  residu_i <- interpolate_sky_points(zenith_dn$sky_points, g,
                                     rmax = ncol(r) / 7)
  sky <- sky_cie - residu_i
  terra::cover(sky, sky_cie)
}
