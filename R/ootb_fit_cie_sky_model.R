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
#'   [sky_grid_segmentation()], see the example. More options translates into
#'   more computing time.
#' @param min_spherical_dist Numeric vector. These values will be passed to the
#'   `min_dist` argument of [vicinity_filter()]. More options translate into
#'   more computing time.
#' @param method Character vector. The methods that will be passed to
#'   [stats::optim()]. More options translates into more computing time.
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
#' r <- caim$Blue
#'
#' com <- compute_complementary_gradients(caim)
#' chroma <- max(com$blue_yellow, com$cyan_red)
#' bin <- apply_thr(chroma, thr_isodata(chroma[!is.na(chroma)]))
#' bin <- bin & apply_thr(com$blue_yellow, -0.2)
#'
#' set.seed(7)
#' gs <- list(
#'   sky_grid_segmentation(z, a, 15, first_ring_different = TRUE)
#' )
#'
#' sky <- ootb_fit_cie_sky_model(r, z, a, m, bin , gs,
#'                               method = "BFGS",
#'                               min_spherical_dist = 0)
#'
#' # set.seed(7)
#' # gs <- list(
#' #   sky_grid_segmentation(z, a, 1.875, first_ring_different = TRUE),
#' #   sky_grid_segmentation(z, a, 6, first_ring_different = TRUE),
#' #   sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
#' # )
#' #
#' # sky <- ootb_fit_cie_sky_model(r, z, a, m, bin , gs,
#' #                               method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
#' #                               min_spherical_dist = seq(0, 9, 3))
#'
#' sky$sky
#' plot(r/sky$sky)
#' sky$model_validation$rmse
#' plot(select_sky_vault_region(r/sky$sky, 0.95, 1.05))
#' plot(sky$model_validation$pred, sky$model_validation$obs)
#' abline(0,1)
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#' }
ootb_fit_cie_sky_model <- function(r, z, a, m, bin, gs,
                                   min_spherical_dist,
                                   method = c("Nelder-Mead", "BFGS",
                                              "CG", "SANN", "Brent")
                                   ) {

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

  sun_zenith_azimuth2 <- extract_sun_zenith_azimuth(r, z, a,
                                              select_sky_vault_region(z, 0, 80),
                                              NULL)

  if (is.null(bin)) {
    # Binarize
    skip_point_extraction_n_filtering <- TRUE

    ## sun mask
    m_sun <- select_circumsolar_region(z, a, sun_zenith_azimuth2, 10)

    ## other mask
    params <- polar_qtree_strict(r, z, a, m, angle_width = 45, diagnose = TRUE)
    seg <- polar_qtree_strict(r, z, a, m,
                              scale_parameter = quantile(params$delta, 0.75),
                              angle_width = 45,
                              min_size_px = mean(params$n_pixels)/4^3,
                              parallel = TRUE,
                              diagnose = FALSE)
    cell_size <- extract_feature(seg, seg, length)
    m2 <- cell_size > mean(params$n_pixels)/4^2

    ## calc thr
    thrs <- calc_thrs(r, z, a,  m & !m_sun & m2,
                      angle_width = 30,
                      fov = 60,
                      method = "thr_isodata")
    thrs_f <- filter_thrs(thrs, r, z, a,
                          g = sky_grid_segmentation(z, a, 10),
                          laxity = 2.5,
                          per_direction = TRUE,
                          n_min = 3)
    thr <- interpolate_planar(thrs_f$rr$sky_points, r,
                              rmax = ncol(r)/7, col_id = "dn")
    thr[is.na(thr)] <- mean(thrs_f$rr$sky_points$dn)
    bin <- apply_thr(r, thr)

    # Detect sky dn
    thrs <- calc_thrs(r, z, a,  m & !m_sun & m2,
                      angle_width = 10,
                      fov = 45,
                      method = "thr_twocorner_uc")

  } else {
    skip_point_extraction_n_filtering <- FALSE
  }

  .fun <- function(g) {
    # Sun coordinate
    sun_zenith_azimuth <- extract_sun_zenith_azimuth(r, z, a, bin, g,
                                                     chi_max_sun = 30)

    # Obtain parameter
    dist_to_black <- optim_dist_to_black(r, z, a, m, bin, g)

    if (skip_point_extraction_n_filtering) {
      thrs_f <- filter_thrs(thrs, r, z, a,
                            g = g,
                            laxity = 2.5,
                            per_direction = TRUE,
                            n_min = 0)
      rr <- thrs_f$rr
      i <- sample(1:nrow(rr$sky_points), size = 100)
      rr$sky_points <- rr$sky_points[i,]

      min_spherical_dist <- 0 #overwrite user input

    } else {
      # Sky points
      sky_points <- extract_sky_points(r, bin, g, dist_to_black)

      ## Apply filters
      sky_points <- vicinity_filter(sky_points, r, min_dist = 3)

      cv <- terra::focal(r, 3, sd) / terra::focal(r, 3, mean)
      sky_points <- sor_filter(sky_points, cv, z, a,
                               k = 10,
                               chi_max = 30,
                               laxity = 2,
                               cutoff_side = "right")

      sky_points <- tryCatch(sor_filter(sky_points, r, z, a,
                                        k = 20,
                                        chi_max = 45,
                                        laxity = 3,
                                        cutoff_side = "both",
                                        trend = 3),
                             error = function(e) sor_filter(sky_points, r, z, a,
                                                            k = 5,
                                                            chi_max = 20,
                                                            laxity = 2,
                                                            cutoff_side = "left")
      )

      rr <- extract_rel_radiance(r, z, a, sky_points, no_of_points = 3,
                                 use_window = !is.null(dist_to_black))
    }


    # Model
    optim_methods <- method

    models <- c(
      Map(function(x) fit_cie_sky_model(rr, sun_zenith_azimuth,
                                        twilight = 0,
                                        method = x), optim_methods),
      Map(function(x) fit_cie_sky_model(rr, sun_zenith_azimuth2,
                                        twilight = 90,
                                        method = x), optim_methods)
    )

    metric <- Map(.get_metric, models)
    i <- which.min(metric)
    model <- models[[i]]
    method <- c(optim_methods, optim_methods)[i]


    # Optim sun
    delta_chi <- 1e10
    while (delta_chi > pi/180) {
      pre_optim_sun <- model$sun_zenith_azimuth
      sun <- Map(function(x) optim_sun_zenith_azimuth(model$sun_zenith_azimuth,
                                                      rr,
                                                      model$coef,
                                                      method = x), optim_methods)

      u <- Map(function(x) !is.na(x[1]), sun) %>% unlist()
      if (any(u)) {
        sun <- Map(function(i) sun[[i]], seq_along(optim_methods)[u])
        models <- Map(function(x) {
          fit_cie_sky_model(rr, unname(x),
                            custom_sky_coef = model$coef,
                            twilight = 90,
                            method = method)
        }, sun)
        metric <- Map(.get_metric, models)
        i <- which.min(metric)
        model <- models[[i]]
      }
      delta_chi <- calc_spherical_distance(pre_optim_sun[1],
                                           pre_optim_sun[2],
                                           model$sun_zenith_azimuth[1],
                                           model$sun_zenith_azimuth[2])
    }
    sun_zenith_azimuth <- model$sun_zenith_azimuth

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

    if (skip_point_extraction_n_filtering) {
      sky_points <- extract_sky_points(r, bin, g, dist_to_black)
    } else {
      sky_points <- rr$sky_points[, c("row", "col")]
    }

    list(sky = .get_sky_cie(z, a, model),
         model = model,
         model_validation = model_validation,
         dist_to_black = dist_to_black,
         use_window = !is.null(dist_to_black),
         min_spherical_dist = min_spherical_dist,
         sky_points = sky_points,
         bin = bin
         )

  }


  # Try grids
  skies <- Map(function(g) {
    tryCatch(.fun(g), error = function(e) NULL)
  }, gs)

  metric <- Map(function(x) tryCatch(.get_metric(x$model_validation),
                                     error = function(e) 1e10
                                     ), skies) %>% unlist()
  i <- which.min(metric)
  if (metric[i] == 1e10) {
    warning("`ootb_fit_cie_sky_model()` failed unexpectedly.")
    sky <- NULL
  } else {
    sky <- skies[[i]]
    sky$g <- gs[[i]]
  }
  sky
}
