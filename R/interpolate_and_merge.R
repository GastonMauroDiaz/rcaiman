#' Interpolate sky data into a raster and merge it with a sky model raster
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. See [ootb_sky_reconstruction()] for
#' further details.
#'
#' @inheritParams ootb_mblt
#' @inheritParams extract_rl
#' @param ootb_sky An object of the class `list` that is the result of calling
#'   [ootb_sky_reconstruction()].
#' @param expand Logical vector of length one. If `TRUE`, [expand_sky_points()]
#'   will be applied internally.
#'
#' @family Sky Reconstruction Functions
#'
#' @return An object of class [SpatRaster-class].
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' bsi <- (caim$Blue-caim$Red) / (caim$Blue+caim$Red)
#'
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r <- caim$Blue
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' sun_coord$zenith_azimuth
#'
#' .a <- azimuth_image(z, orientation = sun_coord$zenith_azimuth[2]+90)
#' seg <- sectors_segmentation(.a, 180) * rings_segmentation(z, 30)
#' bin <- regional_thresholding(r, seg, method = "thr_isodata")
#'
#' mx <- optim_normalize(caim, bin)
#' caim <- normalize(caim, mx = mx, force_range = TRUE)
#' ecaim <- enhance_caim(caim, m, HSV(239, 0.85, 0.5))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m])) & bsi > 0.1
#'
#' set.seed(7)
#' g <- sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
#' sky <- ootb_sky_reconstruction(r, z, a, m, bin & mask_hs(z, 0, 80), g,
#'                                sor_filter_cv = TRUE, sor_filter_dn = TRUE,
#'                                min_spherical_dist = 3)
#'
#' ratio <- r/sky$sky
#' .r <- normalize(ecaim *
#'                   normalize(ratio, 0, 1, TRUE) * normalize(bsi, 0, 1, TRUE))
#' bin <- apply_thr(.r, thr_isodata(.r[m]))
#'
#' g <- sky_grid_segmentation(z, a, 5, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_plant = 3,
#'                                  min_raster_dist = 9)
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 4, pch = 10)
#'
#' plot(sky$sky)
#' sky2 <- interpolate_and_merge(r, z, a, sky_points, sky)
#' plot(sky2)
#' calc_oor_index(r, sky2)
#' sky$oor_index
#' }
interpolate_and_merge <- function(r, z, a, sky_points, ootb_sky,
                                  expand = TRUE) {
  reg <- ootb_sky$model_validation$reg
  model <- ootb_sky$model
  error <- stats::median(abs(1 - reg$model$y/reg$model$x))
  w <-  1 - stats::plogis(error, 0.042, 0.05) + stats::plogis(0, 0.042, 0.05)
  r2 <- summary(reg) %>% .$r.squared
  if (!is.numeric(r2)) r2 <- 0
  k <- 2 + round(r2 * 10)
  p <- abs(model$coef[1]) + log(model$coef[3]+1)
  p <- 2 + sqrt(p)

  i <- grDevices::chull(sky_points$col, nrow(r) - sky_points$row)
  obs <- extract_dn(r, sky_points[i, ])[,3]
  est <- extract_dn(ootb_sky$sky, sky_points[i, ])[,3]
  res <- obs - est

  rmax <- pmax(ncol(r)/14,
               ncol(r)/7 * (1 - stats::median(res)/max(obs))) %>%
    round()

  sky_cie <- .get_sky_cie(z, a, model)

  if (expand) {
    sky_points <- expand_sky_points(r, z, a, sky_points, k = k, p = p)
    rmax <- rmax/2
    sky <- interpolate_sky_points(sky_points, r, k = k, p = p,
                                  rmax = rmax, col_id = "dn")
  } else {
    sky_points <- extract_dn(r, sky_points)
    sky <- interpolate_sky_points(sky_points, r, k = k, p = p,
                                  rmax = rmax, col_id = 3)
  }


  sky <- sky * (1 - w) + sky_cie * w
  sky <- terra::cover(sky, sky_cie)
  names(sky) <- paste0("Weighted average, ",
                       "w=", round(w, 2), ", ",
                       "k=", k, ", ",
                       "p=", round(p, 2), ", ",
                       "rmax=", rmax)

  sky
}


# param weighted_average A logical vector of length one. If `FALSE`, the
#   residuals of the CIE sky model (\eqn{residuals=model-data}) are
#   interpolated instead of the sky digital numbers (the data), and the final
#   sky reconstruction is obtained by subtracting the interpolated residuals to
#   the CIE sky model instead of by calculating a weighted average
#   parameterized by the user, as when the method proposed by
#   \insertCite{Lang2010;textual}{rcaiman} are followed.
# residu <- sky_cie - r
# sky_points <- extract_dn(residu, sky_points)
# residu_i <- interpolate_sky_points(sky_points, r, k = k, p = p,
#                                    rmax = rmax, col_id = 3)
# sky <- sky_cie - residu_i
# sky <- terra::cover(sky, sky_cie)
# names(sky) <- paste0("Residual interpolation, ",
#                      "k=", k, ", ",
#                      "p=", round(p, 2), ", ",
#                      "rmax=", rmax/3)
