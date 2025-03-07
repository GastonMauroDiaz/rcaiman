#' Interpolate sky data into a raster and merge it with a sky model raster
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly
#' presenting and testing this pipeline is under preparation.
#'
#' @inheritParams ootb_mblt
#' @param sky_points An object of class `data.frame`. The output of
#'   [extract_sky_points()] or [expand_sky_points].
#' @param ootb_sky An object of the class `list` that is the result of calling
#'   [ootb_fit_cie_sky_model()].
#' @param rmax_tune Numeric vector of length one. It must be a positive integer.
#'   It is used to fine tune the `rmax` argument that is computed internally
#'   (see Details).
#'
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
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path))
#'
#' sky <- interpolate_and_merge(r, z, a, ootb_sky$sky_points[,c("row", "col")],
#'                              ootb_sky)
#'
#' bin <- apply_thr(r/sky$sky, 0.75)
#' plot(bin)
#'
#' g <- sky_grid_segmentation(z, a, 5, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_plant = 3,
#'                                  min_raster_dist = 9)
#'
#' points(ootb_sky$sky_points$col, nrow(caim) - ootb_sky$sky_points$row,
#'        col = 2, pch = 10)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 4, pch = 10)
#' sky_points <- extract_rl(r, z, a, sky_points, no_of_points = NULL)$sky_points
#' i <- sor_filter(sky_points, r, k = 3)
#'
#' sky2 <- interpolate_and_merge(r, z, a, sky_points[i, c("row", "col")],
#'                               ootb_sky, rmax_tune = 0.75)
#' plot(sky2$sky)
#'
#' sky_points <- expand_sky_points(r, z, a, sky_points[i, c("row", "col")],
#'                                 angle_width = 1.875,
#'                                 k= 3, rmax = 20)
#' sky_points <- sor_filter(sky_points) %>% sky_points[. ,]
#' sky2 <- interpolate_and_merge(r, z, a, sky_points, ootb_sky, rmax_tune = 0.75)
#' plot(sky2$sky)
#' }
interpolate_and_merge <- function(r, z, a, sky_points, ootb_sky,
                                  rmax_tune = 1) {
  model <- ootb_sky$model
  sky_cie <- cie_sky_model_raster(z, a,
                                  model$sun_coord$zenith_azimuth,
                                  model$coef) * model$zenith_dn
  names(sky_cie) <- "CIE sky"

  r2 <- ootb_sky$model_validation$r_squared
  error <- stats::median(abs(1 - ootb_sky$model_validation$observed /
                               ootb_sky$model_validation$predicted))

  expand <- FALSE
  if (!is.null(sky_points$initial)) {
    expand <- TRUE
    sky_points2 <- sky_points
    sky_points <- sky_points[sky_points$initial, c("row", "col")]
  }

  w <-  1 - stats::plogis(error, 0.042, 0.05) + stats::plogis(0, 0.042, 0.05)
  if (!is.numeric(r2)) r2 <- 0
  k <- 2 + round(r2 * 10)
  p <- abs(model$coef[1]) + log(model$coef[3]+1)
  p <- 2 + sqrt(p)

  i <- grDevices::chull(sky_points$col, nrow(r) - sky_points$row)
  obs <- extract_dn(r, sky_points[i, ])[,3]
  est <- extract_dn(sky_cie, sky_points[i, ])[,3]
  res <- obs - est

  rmax <- pmax(ncol(r)/14,
               ncol(r)/7 * (1 - stats::median(res)/max(obs))) %>%
    round()

  if (expand) {
    sky_points <- sky_points2
    sky <- interpolate_sky_points(sky_points, r, k = k, p = p,
                                  rmax = rmax * rmax_tune, col_id = "dn")
  } else {
    sky_points <- extract_dn(r, sky_points)
    sky <- interpolate_sky_points(sky_points, r, k = k, p = p,
                                  rmax = rmax * rmax_tune, col_id = 3)
  }

  sky <- sky * (1 - w) + sky_cie * w
  sky <- terra::cover(sky, sky_cie)
  names(sky) <- paste0("Weighted average, ",
                       "w=", round(w, 2), ", ",
                       "k=", k, ", ",
                       "p=", round(p, 2), ", ",
                       "rmax=", rmax)

  list(sky = sky,
       w = w,
       k=  k,
       p = p,
       rmax = rmax * rmax_tune)
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
