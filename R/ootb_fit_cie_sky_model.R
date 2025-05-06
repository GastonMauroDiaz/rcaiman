#' Out-of-the-box fit CIE sky model
#'
#' Out-of-the-box fit CIE sky model
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly presenting and
#' testing this pipeline is under preparation.
#'
#' @inheritParams sky_grid_segmentation
#' @inheritParams extract_sky_points
#' @inheritParams fit_trend_surface
#' @inheritParams fit_cie_sky_model
#' @param m [SpatRaster-class]. A mask, check [select_sky_vault_region()].
#' @param gs An object of the class _list_. A list with the output of
#'   [sky_grid_segmentation()], see the example.
#' @param min_spherical_dist Numeric vector. These values will be passed to the
#'   `min_dist` argument of [vicinity_filter()].
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @return A _list_ with the following components:
#' \itemize{
#'   \item An object of the [SpatRaster-class]) with the predicted digital
#'   number values.
#'   \item The output produced by [fit_cie_sky_model()].
#'   \item The output produced by [validate_cie_sky_model()].
#'   \item The `dist_to_black` argument used in [extract_sky_points()].
#'   \item The `use_window` argument used in [extract_rel_radiance()].
#'   \item The `min_spherical_dist` argument used to filter the sky points.
#'   \item The `sky_points` argument used in [extract_rel_radiance()].
#' }
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r <- caim$Blue
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_zenith_azimuth <- extract_sun_zenith_azimuth(r, z, a, bin, g)
#'
#' .a <- azimuth_image(z, orientation = sun_zenith_azimuth[2] + 90)
#' seg <- sectors_segmentation(.a, 180) * rings_segmentation(z, 30)
#' bin <- regional_thresholding(r, seg, method = "thr_isodata")
#'
#' mx <- optim_max(caim, bin)
#' caim <- normalize_minmax(caim, mx = mx, force_range = TRUE)
#' ecaim <- enhance_caim(caim, m, polarLAB(50, 17, 293))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' bin <- bin & select_sky_vault_region(z, 0, 85)
#' plot(bin)
#'
#' set.seed(7)
#' gs <- list(
#'   sky_grid_segmentation(z, a, 1.875, first_ring_different = TRUE),
#'   sky_grid_segmentation(z, a, 6, first_ring_different = TRUE),
#'   sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
#' )
#'
#' sky <- ootb_fit_cie_sky_model(r, z, a, m, bin , gs,
#'                               min_spherical_dist = seq(3, 9, 3))
#'
#' sky$sky
#' plot(sky$sky)
#' sky$model_validation$rmse
#' plot(select_sky_vault_region(r/sky$sky, 0.95, 1.05))
#' plot(sky$model_validation$pred, sky$model_validation$obs)
#' abline(0,1)
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#' }
ootb_fit_cie_sky_model <- function(r, z, a, m, bin, gs,
                                   min_spherical_dist) {

  .is_single_layer_raster(m, "m")
  .is_logic_and_NA_free(m)
  stopifnot(compareGeom(r, m) == TRUE)

  .get_metric <- function(model) {
    #10.2134/agronj2003.1442
    x <- model$pred
    y <- model$obs
    reg <- lm(x~y)
    m <- stats::coef(reg)[2]
    r_squared <- summary(reg) %>% .$r.squared
    SB <- (mean(x) - mean(y))^2
    NU <- (1 - m)^2 * mean(x^2)
    LC <- (1 - r_squared) * mean(y^2)
    MSE <- SB + NU + LC
    MSE
  }

  .get_sky_cie <- function(z, a, model) {
    sky_cie <- cie_sky_image(z, a,
                             model$sun_zenith_azimuth,
                             model$coef) * model$zenith_dn
    names(sky_cie) <- "CIE sky"
    sky_cie
  }

  .fun <- function(g) {
    # Sun coordinate
    sun_zenith_azimuth <- extract_sun_zenith_azimuth(r, z, a, bin, g,
                                                     max_angular_dist = 30,
                                                     method = "object-based")
    sun_zenith_azimuth2 <- extract_sun_zenith_azimuth(r, z, a, bin, g,
                                                      max_angular_dist = 30,
                                                      method = "pixel-based")

    # Sky points
    dist_to_black <- optim_dist_to_black(r, z, a, m, bin, g)
    sky_points <- extract_sky_points(r, bin, g, dist_to_black)

    ## Apply filters
    sky_points <- vicinity_filter(sky_points, r, min_dist = 3)

    cv <- terra::focal(r, 3, sd) / terra::focal(r, 3, mean)
    sky_points <- sor_filter(sky_points, cv, z, a,
                             k = 10,
                             rmax = 30,
                             thr = 2,
                             cutoff_side = "right")

    sky_points <- tryCatch(sor_filter(sky_points, r, z, a,
                                      k = 20,
                                      rmax = 45,
                                      thr = 3,
                                      cutoff_side = "both",
                                      trend = 3),
                           error = function(e) sor_filter(sky_points, r, z, a,
                                                          k = 5,
                                                          rmax = 20,
                                                          thr = 2,
                                                          cutoff_side = "left")
    )


    # Model
    rr <- extract_rel_radiance(r, z, a, sky_points, no_of_points = 3,
                               use_window = !is.null(dist_to_black))

    methods <- c("Nelder-Mead", "BFGS", "CG", "SANN", "Brent")

    models <- c(
              Map(function(x) fit_cie_sky_model(rr, sun_zenith_azimuth,
                                                twilight = 0,
                                                method = x), methods),
              Map(function(x) fit_cie_sky_model(rr, sun_zenith_azimuth2,
                                                twilight = 90,
                                                method = x), methods)
              )

    metric <- Map(.get_metric, models)
    i <- which.min(metric)
    model <- models[[i]]
    method <- methods[i]
    sun_zenith_azimuth <- model$sun_zenith_azimuth


    # Optim sun
    sun <- Map(function(x) optim_sun_zenith_azimuth(sun_zenith_azimuth,
                                                    rr,
                                                    model$coef,
                                                    method = x), methods)

    u <- Map(function(x) !is.na(x[1]), sun) %>% unlist()
    if (any(u)) {
      sun <- Map(function(i) sun[[i]], seq_along(methods)[u])
      models <- Map(function(x) {
        fit_cie_sky_model(rr, x,
                          custom_sky_coef = model$coef,
                          twilight = 90,
                          method = method)
      }, sun)
      metric <- Map(.get_metric, models)
      i <- which.min(metric)
      model <- models[[i]]
    }

    # Try subsamples
    if (any(min_spherical_dist != 0)) {
      .get_rrs <- function(min_spherical_dist) {
        sky_points <- vicinity_filter(sky_points, r, z, a, min_spherical_dist,
                                      use_window = !is.null(dist_to_black))
        extract_rel_radiance(r, z, a, sky_points, no_of_points = 3,
                                      use_window = !is.null(dist_to_black))
      }
      rrs <- Map(.get_rrs, min_spherical_dist)
      models <- Map(function(rr) {
        fit_cie_sky_model(rr, sun_zenith_azimuth,
                          custom_sky_coef = model$coef,
                          twilight = 90,
                          method = method)
      }, rrs)
      metric <- Map(.get_metric, models)
      i <- which.min(metric)
      model <- models[[i]]
      rr <- rrs[[i]]
      min_spherical_dist <- min_spherical_dist[i]
    } else {
      model <- fit_cie_sky_model(rr, sun_zenith_azimuth,
                                 custom_sky_coef = model$coef,
                                 twilight = 90,
                                 method = method)
    }


    model_validation <- validate_cie_sky_model(model, rr, k = 10)

    list(sky = .get_sky_cie(z, a, model),
         model = model,
         model_validation = model_validation,
         dist_to_black = dist_to_black,
         use_window = !is.null(dist_to_black),
         min_spherical_dist = min_spherical_dist,
         sky_points = rr$sky_points[, c("row", "col")]
         )

  }


  # Try grids
  skies <- Map(function(g) {
    tryCatch(.fun(g),
             error = function(e) list(model_validation = list(pred = 0,
                                                              obs = 1e10))
             )
  }, gs)

  metric <- Map(function(x) .get_metric(x$model_validation), skies) %>% unlist()
  i <- which.min(metric)
  if (metric[i] == 1e10) {
    sky <- NULL
  } else {
    sky <- skies[[i]]
    sky$g <- gs[[i]]
  }
  sky
}
