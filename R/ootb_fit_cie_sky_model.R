#' Out-of-the-box fit CIE sky model
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions [extract_sun_coord()], [extract_sky_points()], [sor_filter()],
#' [extract_rl()], [fit_cie_sky_model()], and [validate_cie_sky_model()].
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly presenting and
#' testing this pipeline is under preparation.
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_trend_surface
#' @inheritParams extract_sky_points
#' @param m [SpatRaster-class]. A mask, check [mask_hs()].
#' @param refine_sun_coord Logical vector of length one
#' @param sor_filter_cv Logical vector of length one. If `TRUE`, [sor_filter()]
#'   will be applied internally to filter out points using local variability at
#'   the scale of the \eqn{3 \times 3} window using for data extraction as the
#'   criterion. Local means not farther than 30 degrees.
#' @param sor_filter_dn Logical vector of length one. If `TRUE`, [sor_filter()]
#'   will be applied internally to filter out points using local darkness as the
#'   criterion. Local means not farther than 20 degrees.
#' @param input_sky_points An object of class *data.frame* with the same
#'   structure than the output of [extract_sky_points()]. This argument is
#'   convinient to provide manually digitized points, see [fit_cie_sky_model()]
#'   for details.
#' @inheritParams extract_rl
#'
#' @export
#'
#' @family Sky Reconstruction Functions
#'
#' @references \insertAllCited{}
#'
#' @return A _list_ with the following components:
#' \itemize{
#'   \item An object of the [SpatRaster-class]) with the predicted digital
#'   number values.
#'   \item The output produced by [fit_cie_sky_model()].
#'   \item The output produced by [validate_cie_sky_model()].
#'   \item The `dist_to_black` argument used in [fit_cie_sky_model()].
#'   \item The `sky_points` argument used in [extract_rl()].
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
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' sun_coord$zenith_azimuth
#'
#' .a <- azimuth_image(z, orientation = sun_coord$zenith_azimuth[2]+90)
#' seg <- sectors_segmentation(.a, 180) * rings_segmentation(z, 30)
#' bin <- regional_thresholding(r, seg, method = "thr_isodata")
#'
#' mx <- optim_normalize(caim, bin)
#' caim <- normalize(caim, mx = mx, force_range = TRUE)
#' ecaim <- enhance_caim(caim, m, polarLAB(50, 17, 293))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' bin <- bin & mask_hs(z, 0, 85)
#' plot(bin)
#'
#' set.seed(7)
#' g <- sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
#' sky <- ootb_fit_cie_sky_model(r, z, a, m, bin , g,
#'                               sor_filter_cv = TRUE, sor_filter_dn = TRUE,
#'                               refine_sun_coord = TRUE,
#'                               min_spherical_dist = 3)
#'
#' sky$sky
#' plot(sky$sky)
#' sky$model_validation$rmse
#' plot(r/sky$sky>1.15)
#' plot(sky$model_validation$predicted, sky$model_validation$observed)
#' abline(0,1)
#' error <- sky$model_validation$predicted - sky$model_validation$observed
#' plot(sky$sky_points$z[!sky$sky_points$is_outlier], error,
#'      xlab = "zenith angle", ylab = "relative radiance error")
#' abline(h = 0)
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#'
#' }
ootb_fit_cie_sky_model <- function(r, z, a, m, bin, g,
                                   sor_filter_cv = FALSE,
                                   sor_filter_dn = FALSE,
                                   refine_sun_coord = FALSE,
                                   min_spherical_dist = NULL,
                                   input_sky_points = NULL
                                   ) {

  .is_single_layer_raster(m, "m")
  .is_logic_and_NA_free(m)
  stopifnot(compareGeom(r, m) == TRUE)

  .get_metric <- function(model) {
    stats::median(abs(model$pred - model$obs))
  }

  dist_to_black <- optim_dist_to_black(r, z, a, m, bin, g)

  sky_points <- extract_sky_points(r, bin, g, dist_to_black)
  if (!is.null(input_sky_points)) {
    sky_points <- rbind(input_sky_points, sky_points)
  }

  sun_coord <- extract_sun_coord(r, z, a, bin, g)

  rl <- extract_rl(r, z, a, sky_points, use_window = dist_to_black != 1,
                   min_spherical_dist = min_spherical_dist)
  if (sd(rl$sky_points$rl) < 0.01) {
    warning("Overexposed image")
  }

  methods <- c("Nelder-Mead", "BFGS", "CG", "SANN", "Brent")
  models <- Map(function(x) fit_cie_sky_model(rl, sun_coord,
                                              twilight = 60,
                                              method = x), methods)
  metric <- Map(.get_metric, models)
  i <- which.min(metric)
  model <- models[[i]]
  method <- methods[i]
  sun_coord <- model$sun_coord

  # START sampling on the out-of-range zone ####
  if (dist_to_black != 1) {
    model.2 <- model
    model$pred <- model$obs * 1e+10
    rl.2 <- rl
    while (.get_metric(model) - .get_metric(model.2) > 0.0001) {
      model <- model.2
      rl <- rl.2

      ratio <- r / .get_sky_cie(z, a, model)
      ratio[is.infinite(ratio)] <- 1e+10
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
      out.of.range_ratio[is.na(out.of.range_ratio)] <- 0

      v <- terra::cellFromRowCol(r, rl$sky_points$row, rl$sky_points$col) %>%
        xyFromCell(r, .) %>% vect()
      v <- terra::buffer(v, dist_to_black)
      bin[v] <- 0
      bin <- bin & out.of.range_ratio != 0
      if (calc_co(bin, z, a, m) > 0.05) {
        sky_points.2 <- extract_sky_points(out.of.range_ratio, bin, g,
                                           dist_to_black + 2)
        sky_points.2 <- rbind(rl$sky_points[ , c("row", "col")], sky_points.2)
        rl.2 <- extract_rl(r, z, a, sky_points.2, g, no_of_points = NULL,
                           use_window = dist_to_black != 1,
                           min_spherical_dist = min_spherical_dist)
        rl.2$sky_points$rl <- rl.2$sky_points$rl/rl$zenith_dn
        rl.2$zenith_dn <- rl$zenith_dn

        model.2 <- fit_cie_sky_model(rl.2, sun_coord,
                                     custom_sky_coef = model$coef,
                                     twilight = 90,
                                     method = method)
      }
    }
  }
  # END sampling on the out-of-range zone ####

  # START thin points ####
  .filter <- function(ds, col_names, thr) { #from extract_sky_points()
    d <- as.matrix(stats::dist(ds[, col_names]))
    indices <- c()
    i <- 0
    while (i < nrow(d)) {
      i <- i + 1
      indices <- c(indices, row.names(d)[i]) #include the point itself (p)
      x <- names(d[i, d[i,] <= thr])
      if (!is.null(x)) {
        # this exclude from future search all the points near p,
        # including itself
        rows2crop <- (1:nrow(d))[match(x, rownames(d))]
        cols2crop <- (1:ncol(d))[match(x, colnames(d))]
        d <- d[-rows2crop, -cols2crop]
      }
      if (is.vector(d)) d <- matrix(d)
    }
    ds[indices,]
  }
  rl$sky_points <- .filter(rl$sky_points, c("col", "row"), 0.1) #rem identicals

  if (sor_filter_cv == TRUE) {
    cv <- terra::focal(r, 3, sd) / terra::focal(r, 3, mean)
    u <- sor_filter(rl$sky_points, cv,
                    k = 10,
                    rmax = 30,
                    thr = 2,
                    cutoff_side = "right")
  } else {
    u <- rep(TRUE, nrow(sky_points))
  }
  if (sor_filter_dn == TRUE) {
    v <- sor_filter(rl$sky_points,
                    k = 5,
                    rmax = 20,
                    thr = 2,
                    cutoff_side = "left")
  } else {
    v <- rep(TRUE, nrow(sky_points))
  }

  i <- u & v
  if (!is.null(input_sky_points)) {
    i[1:nrow(input_sky_points)] <- TRUE
  }
  rl$sky_points <- rl$sky_points[i, ]

  model <- fit_cie_sky_model(rl, sun_coord,
                             custom_sky_coef = model$coef,
                             twilight = 90,
                             method = method)
  # END thin points ####

  # START sun_coord refinement ####
  if (refine_sun_coord) {
    .refine_sun_coord <- function(param) {
      zenith <- param[1] * 9
      azimuth <- param[2] * 36
      pred <- .cie_sky_model(AzP = rl$sky_points$a %>% .degree2radian(),
                             Zp = rl$sky_points$z %>% .degree2radian(),
                             AzS = azimuth %>% .degree2radian(),
                             Zs =  zenith %>% .degree2radian(),
                             model$coef[1], model$coef[2], model$coef[3],
                             model$coef[4], model$coef[5])
      .get_metric(list(pred = pred, obs = model$obs))
    }
    try(fit <- stats::optim(c(sun_coord$zenith_azimuth[1]/9,
                              sun_coord$zenith_azimuth[2]/36),
                            .refine_sun_coord,
                            lower = 0,
                            upper = 10,
                            method = "L-BFGS-B"), silent = TRUE)
    try(sun_coord$zenith_azimuth <- c(fit$par[1]*9, fit$par[2]*36),
        silent = TRUE)
    model <- fit_cie_sky_model(rl, sun_coord,
                               custom_sky_coef = model$coef,
                               twilight = 90,
                               method = method)
  }
  # END sun coord refinement ####

  models <- Map(function(x) fit_cie_sky_model(rl, sun_coord, method = x,
                                              twilight = 90),
                methods)
  metric <- Map(.get_metric, models)
  i <- which.min(metric)
  model <- models[[i]]

  model_validation <- validate_cie_sky_model(r, z, a, rl, model,
                                             use_window = dist_to_black != 1,
                                             k = 10)

  rl$sky_points$is_outlier <- model_validation$is_outlier

  list(sky = .get_sky_cie(z, a, model),
       model = model,
       model_validation = model_validation,
       dist_to_black = dist_to_black,
       sky_points = rl$sky_points
       )
}
