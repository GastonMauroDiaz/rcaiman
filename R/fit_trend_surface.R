#' Fit a trend surface to sky digital numbers
#'
#' Fit a trend surface using [spatial::surf.ls()] as workhorse function.
#'
#' This function is meant to be used after [fit_coneshaped_model()].
#'
#' The first application of trend surface fitting to modelling sky digital
#' numbers was presented in \insertCite{Diaz2018;textual}{rcaiman}, under the
#' heading *Estimation of the sky DN as a previous step for our method*. The
#' example shows a pipeline that resemble the one presented in that paper. The
#' original idea was developed after the programming of [extract_sky_points()].
#'
#'
#' @inheritParams interpolate_planar
#' @inheritParams spatial::surf.ls
#'
#' @note If an incomplete above-canopy image is available as filling source,
#'   non-sky pixels should be turned `NA` or they will be erroneously considered
#'   as sky pixels.
#'
#' @return A list with an object of class [SpatRaster-class] and of class `trls`
#'   (see [spatial::surf.ls()]).
#' @export
#'
#' @seealso [thr_mblt()]
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' mx <- optim_max(caim, bin)
#' caim <- normalize_minmax(caim, 0, mx, TRUE)
#'
#' sky_blue <- polarLAB(50, 17, 293)
#' ecaim <- enhance_caim(caim, m, sky_blue = sky_blue)
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#'
#' g <- sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_black = 3)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' rr <- extract_rel_radiance(r, z, a, sky_points)
#'
#' sky_s <- fit_trend_surface(rr$sky_points, r, col_id = "dn")
#' plot(sky_s$raster)
#' plot(r/sky_s$raster)
#' apply_thr(r/sky_s$raster, 0.5) %>% plot()
#'
#' model <- fit_coneshaped_model(rr$sky_points)
#' summary(model$model)
#' sky_cs <- model$fun(z, a) * rr$zenith_dn
#' plot(sky_cs)
#'
#' .r <- r
#' .r[!bin] <- sky_cs
#' plot(.r)
#' sky_points <- extract_sky_points(.r, m, g, dist_to_black = 3)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' rr <- extract_rel_radiance(.r, z, a, sky_points)
#'
#' sky_s <- fit_trend_surface(rr$sky_points, .r, col_id = "dn")
#' plot(sky_s$raster)
#' plot(r/sky_s$raster)
#' apply_thr(r/sky_s$raster, 0.5) %>% plot()
#' }
fit_trend_surface <- function(sky_points,
                              r,
                              np = 6,
                              col_id = "dn") {

  stopifnot(length(np) == 1)
  xy <- terra::cellFromRowCol(r, sky_points$row, sky_points$col) %>%
    terra::xyFromCell(r, .)
  fit <- spatial::surf.ls(x = xy[, 1],
                          y = xy[, 2],
                          z = sky_points[, col_id],
                          np)
  xl <- terra::xmin(r)
  xu <- terra::xmax(r)
  yl <- terra::ymin(r)
  yu <- terra::ymax(r)

  out <- spatial::trmat(fit, xl, xu, yl, yu, max(ncol(r), nrow(r)))
  out <- terra::rast(out$z) %>% t %>% flip
  terra::crs(out) <- terra::crs(r)
  terra::ext(out) <- terra::ext(r)
  out <- terra::resample(out, r)

  xy <- xy[grDevices::chull(xy),]
  v <- terra::vect(list(xy), type = "polygon", crs = terra::crs(r))
  v <- terra::rasterize(v, r) %>% is.na()
  out[v] <- NA

  list(raster = out, model = fit)
}
