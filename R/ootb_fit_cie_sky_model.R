#' Out-of-the-box fit CIE sky model
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions [extract_sun_coord()], [extract_sky_points()], [sor_filter()],
#' [extract_rl()] and [fit_cie_sky_model()].
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly presenting and
#' testing this pipeline is under preparation.
#'
#' The validation of the CIE sky model is done with a k-fold approach (k = 10)
#' following \insertCite{Pineiro2008;textual}{rcaiman}
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
#'   \item A _list_ containing the result of model validation:
#'     \itemize{
#'       \item An object of class `lm` (see [stats::lm()]).
#'       \item the coefficient of determination (\eqn{r^2}).
#'       \item predicted values.
#'       \item observed vales.
#'       \item The root mean squared error (RMSE).
#'     }
#'   \item The `dist_to_plant` argument used in [fit_cie_sky_model()].
#'   \item The `sky_points` argument used in [extract_rl()].
#' }
#'
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
#' ecaim <- enhance_caim(caim, m, HSV(239, 0.85, 0.5))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
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
#' plot(sky$sky_points$z[!sky$sky_points$outliers], error,
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

  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  .noise <- function(w = 1) {
    .sd <- apply((skies[, 1:5]), 2, sd) * w
    Map(function(i) stats::rnorm(1, 0, .sd[i]), 1:5) %>% unlist()
  }

  .get_metric <- function(model) {
    stats::median(abs(model$pred - model$obs))
  }

  if (is.null(bin) & (is.null(input_sky_points) | is.logical(refine_sun_coord)))
  {
    stop("When no preliminary binarization is provided through
          the 'bin' argument, then sun coordinates should be provided through
          the 'refine_sun_coord' argument and sky points through the
          'input_sky_points' argunment. See fit_cie_sky_model() for
          more details.")
  }

  if (!is.null(bin)) {
    if (!is.null(g)) {
      # START optimize dist_to_plant ####
      g30 <- sky_grid_segmentation(z, a, 30)
      g30[!m] <- 0
      dist_to_plant <- 11
      sampling_pct <- 0
      while (sampling_pct < 100 & dist_to_plant > 3) {
        dist_to_plant <- dist_to_plant-2
        sky_points <- extract_sky_points(r, bin, g,
                                         dist_to_plant = dist_to_plant)
        v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
          xyFromCell(r, .) %>% vect()
        sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
          (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
      }
      if (sampling_pct < 75) {
        dist_to_plant <- 1
      }
      # END optimize dist_to_plant ####
      if (is.logical(refine_sun_coord)) {
        sun_coord <- extract_sun_coord(r, z, a, bin, g)
      } else {
        sun_coord <- refine_sun_coord
        refine_sun_coord <- FALSE
      }
      sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
    } else {
      if (is.logical(refine_sun_coord)) {
        g10 <- sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
        sun_coord <- extract_sun_coord(r, z, a, bin, g10)
      } else {
        sun_coord <- refine_sun_coord
        refine_sun_coord <- FALSE
      }
      dist_to_plant <- 3
      sky_points <- extract_sky_points(r, bin, g, dist_to_plant = dist_to_plant)
    }
  } else {
    sky_points <- input_sky_points
    input_sky_points <- NULL
    refine_sun_coord <- FALSE
    dist_to_plant <- 3
  }


  if (!is.null(input_sky_points)) {
    sky_points <- rbind(input_sky_points, sky_points)
  }
  rl <- extract_rl(r, z, a, sky_points, use_window = dist_to_plant != 1,
                   min_spherical_dist = min_spherical_dist)
  if (sd(rl$sky_points$rl) < 0.01) {
    warning("Overexposed image")
  }

  method <- c("BFGS", "CG", "SANN")
  models <- Map(function(x) fit_cie_sky_model(rl, sun_coord,
                                              twilight = TRUE,
                                              method = x), method)
  metric <- Map(.get_metric, models)
  i <- which.min(metric)
  model <- models[[i]]
  method <- method[i]
  sun_coord <- model$sun_coord

  # START sampling on the out-of-range zone ####
  if (dist_to_plant != 1 & !is.null(g) & !is.null(bin)) {
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
      v <- terra::buffer(v, dist_to_plant)
      bin[v] <- 0
      bin <- bin & out.of.range_ratio != 0
      if (calc_co(bin, z, a, m) > 0.05) {
        sky_points.2 <- extract_sky_points(out.of.range_ratio, bin, g,
                                           dist_to_plant + 2)
        sky_points.2 <- rbind(rl$sky_points[ , c("row", "col")], sky_points.2)
        rl.2 <- extract_rl(r, z, a, sky_points.2, g, no_of_points = NULL,
                           use_window = dist_to_plant != 1,
                           min_spherical_dist = min_spherical_dist)
        rl.2$sky_points$rl <- rl.2$sky_points$rl/rl$zenith_dn
        rl.2$zenith_dn <- rl$zenith_dn

        model.2 <- fit_cie_sky_model(rl.2, sun_coord,
                                     custom_sky_coef = model$coef,
                                     twilight = FALSE,
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
    u <- sor_filter(rl$sky_points, k = 10,
                    cv, rmax = 30, thr = 0,
                    cutoff_side = "right")
    if (!is.null(input_sky_points)) {
      u[1:nrow(input_sky_points)] <- TRUE
    }
    rl$sky_points <- rl$sky_points[u, ]
  }
  if (sor_filter_dn == TRUE) {
    u <- sor_filter(rl$sky_points, k = 5,
                    rmax = 20, thr = 2, cutoff_side = "left")
    if (!is.null(input_sky_points)) {
      u[1:nrow(input_sky_points)] <- TRUE
    }
    rl$sky_points <- rl$sky_points[u, ]
  }

  model <- fit_cie_sky_model(rl, sun_coord,
                             custom_sky_coef = model$coef,
                             twilight = FALSE,
                             method = method)
  # END thin points ####

  # START sun_coord refinement ####
  if (refine_sun_coord) {
     .refine_sun_coord <- function(zenith, azimuth) {
      sun_coord$zenith_azimuth <- c(zenith, azimuth)
      model <- fit_cie_sky_model(rl, sun_coord,
                                 custom_sky_coef = model$coef,
                                 twilight = FALSE,
                                 method = method)
      .get_metric(model)
    }
    fit <- bbmle::mle2(.refine_sun_coord,
                       list(zenith = sun_coord$zenith_azimuth[1],
                            azimuth = sun_coord$zenith_azimuth[2]),
                       method = "BFGS")
    sun_coord$zenith_azimuth <- fit@coef
    model <- fit_cie_sky_model(rl, sun_coord,
                               custom_sky_coef = model$coef,
                               twilight = FALSE,
                               method = method)
    sun_coord <- model$sun_coord
  }
  # END sun coord refinement ####

  method <- c("BFGS", "CG", "SANN")
  models <- Map(function(x) fit_cie_sky_model(rl, sun_coord, method = x,
                                              twilight = FALSE),
                method)
  metric <- Map(.get_metric, models)
  i <- which.min(metric)
  model <- models[[i]]

  # START K-fold method ####
  ## k=10 based on https://dl.acm.org/doi/10.5555/1643031.1643047
  k <- 10
  folds <- seq_along(rl$sky_points$row)
  folds <- split(folds, 1:k) %>% suppressWarnings()
  x <- c()
  y <- c()
  for (i in 1:k) {
    rl.2 <- rl
    rl.2$sky_points <- rl$sky_points[-folds[[i]],]
    model.2 <- fit_cie_sky_model(rl.2, sun_coord,
                           custom_sky_coef = model$coef + .noise(0.1),
                           twilight = FALSE,
                           method = method)
    x <- c(x, extract_dn(.get_sky_cie(z, a, model.2)/rl$zenith_dn,
                         rl$sky_points[folds[[i]], c("row", "col")],
                         use_window = dist_to_plant != 1)[,3])
    y <- c(y, extract_dn(r/rl$zenith_dn,
                         rl$sky_points[folds[[i]], c("row", "col")],
                         use_window = dist_to_plant != 1)[,3])
  }

  #following Leys2013 10.1016/j.jesp.2013.03.013
  error <- y - x
  u <- abs((error - stats::median(error)) / stats::mad(error)) < 3
  x <- x[u]
  y <- y[u]

  reg <- lm(x~y) #following Pineiro2008 10.1016/j.ecolmodel.2008.05.006
  # END K-fold method ####

  rl$sky_points$outliers <- !u

  list(sky = .get_sky_cie(z, a, model),
       model = model,
       model_validation = list(lm = reg,
                               predicted = reg$model$x,
                               observed = reg$model$y,
                               r_squared = summary(reg) %>% .$r.squared,
                               rmse = .calc_rmse(reg$model$y - reg$model$x)),
       dist_to_plant = dist_to_plant,
       sky_points = rl$sky_points
       )
}
