#' Out-of-the-box sky reconstruction
#'
#' Build an above canopy image from a single below canopy image.
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions, \code{\link{fit_cie_sky_model}} and
#' \code{\link{interpolate_sky_points}}. The code can be easily inspected by
#' calling \code{ootb_sky_reconstruction} --no parenthesis. Advanced users could
#' use that code as a template.
#'
#' This pipeline is based on \insertCite{Lang2010;textual}{rcaiman}. The main
#' differences between the original method by
#' \insertCite{Lang2010;textual}{rcaiman} and the one implemented here are:
#'
#' \itemize{
#'
#' \item it is fully automatic,
#'
#' \item the residuals of the CIE sky model (\eqn{residuals=model-data}) are
#' interpolated instead of the sky digital numbers (the data), and
#'
#' \item the final sky reconstruction is obtained by subtracting the
#' interpolated residuals to the CIE sky model instead of by calculating a
#' weighted average parameterized by the user.
#'
#' }
#'
#' The recommended input for this function is data pre-processed with the HSP
#' software package \insertCite{Lang2013}{rcaiman}. Please, refer to
#' \code{\link{write_sky_points}} for additional details about HSP and refer to
#' \code{\link{fit_cie_sky_model}} and \code{\link{interpolate_sky_points}} to
#' know why the HSP pre-processing is convenient.
#'
#'
#' Providing a \code{filling source} triggers an alternative pipeline in which
#' the sky is fully reconstructed with \code{\link{interpolate_sky_points}}
#' after a dense sampling (\eqn{1 \times 1} degree cells), which is supported by
#' the fact that sky digital numbers will be available for almost every pixel,
#' either from \code{r} gaps or from the filling source --an exception is a
#' filling source with plenty of \code{NA} values.
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_trend_surface
#'
#' @export
#'
#' @family Sky Reconstruction Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' #JPEG file
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' bin <- ootb_mblt(r, z, a)
#' bin <- bin$bin & mask_hs(z, 0, 80)
#' sky <- ootb_sky_reconstruction(r, z, a, bin)
#' plot(sky)
#' ratio <- r/sky
#' plot(ratio)
#' hist(ratio)
#'
#' #preprocessed with HSP
#' path <- system.file("external/DSCN6342.pgm", package = "rcaiman")
#' r <- rast(path)
#' ext(r) <- ext(0,ncol(r),0,nrow(r))
#' crs(r) <- crs(read_caim())
#' r <- normalize(r)
#' z <- zenith_image(ncol(r), lens())
#' a <- azimuth_image(z)
#' bin <- find_sky_pixels(r, z, a)
#' sky <- ootb_sky_reconstruction(r, z, a, bin)
#' sky <- ootb_sky_reconstruction(r, z, a, bin, sky)
#' ratio <- r/sky
#' ratio[is.na(ratio)] <- 0
#' ratio <- normalize(ratio, force_range = TRUE)
#' plot(ratio)
#' g <- sky_grid_segmentation(z, a, 15)
#' plot(defuzzify(ratio, g))
#' }
ootb_sky_reconstruction <- function(r, z, a, bin,
                                    filling_source = NULL) {
  if (is.null(filling_source)) {
    g <- sky_grid_segmentation(z, a, 10)
    sky_points <- extract_sky_points(r, bin, g)
    rl <- extract_rl(r, z, a, sky_points)
    sun_coord <- extract_sun_coord(r, z, a, bin, g)
    model <- fit_cie_sky_model(r, z, a,
                               rl$sky_points,
                               rl$zenith_dn,
                               sun_coord)
    sky_cie <- model$relative_luminance * model$zenith_dn
    sky_cie <- normalize(sky_cie, 0, 1, TRUE)

    residu <- sky_cie - r
    sky_points <- suppressWarnings(extract_rl(residu, z, a, sky_points, NULL))
    residu_i <- interpolate_sky_points(sky_points$sky_points, g,
                                       rmax = ncol(r)/7)

    sky <- sky_cie - residu_i
    sky <- terra::cover(sky, sky_cie)
  } else {
    g <- sky_grid_segmentation(z, a, 10)
    terra::compareGeom(r, filling_source)
    r[!bin] <- filling_source[!bin]
    sky_points <- extract_sky_points(r, !is.na(z),
                                     sky_grid_segmentation(z, a, 1),
                                     dist_to_plant = NULL,
                                     min_raster_dist = NULL)
    sky_points <- extract_rl(r, z, a, sky_points, NULL,
                             use_window = FALSE)$sky_points
    sky_points <- sky_points[!is.na(sky_points$rl),]
    sky <- interpolate_sky_points(sky_points, g)
  }
  sky
}
