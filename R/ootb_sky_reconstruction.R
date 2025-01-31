#' Out-of-the-box sky reconstruction
#'
#' Build an above canopy image from a single below canopy image
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions [fit_cie_sky_model()] and [interpolate_sky_points()].
#'
#' The pipeline is an automatic version of the
#' \insertCite{Lang2010;textual}{rcaiman} method. A paper for thoroughly
#' presenting and testing this pipeline is under preparation.
#'
#' The __out-of-range index__ is calculated as foollow:
#'
#' \eqn{\sum_{i = 1}^{N}(r_i/sky_i)^2},
#'
#' where \eqn{r} is the `r` argument, \eqn{sky} is the
#' raster obtained from the fitted model with [cie_sky_model_raster()] and
#' `zenith_dn`, \eqn{i} is the index that represents the position of a given
#' pixel on the raster grid, and \eqn{N} is the total number of pixels that
#' satisfy: \eqn{r_i/sky_i<0} or \eqn{r_i/sky_i>1}.
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_trend_surface
#' @inheritParams extract_sky_points
#' @param m [SpatRaster-class]. A mask, check [mask_hs()].
#' @param refine_sun_coord Logical vector of length one
#' @param sor_filter_cv Logical vector of length one
#' @param sor_filter_dn Logical vector of length one
#' @param interpolate Logical vector of length one. If `TRUE`,
#'   [interpolate_sky_points()] will be used.
#' @param input_sky_points An object of class *data.frame* with the same
#'   structure than the result of a call to [extract_sky_points()]. The
#'   [ImageJ](https://imagej.net/ij/) software package can be used to manually
#'   digitize points. See [extract_dn()] for details.
#'
#' @export
#'
#' @family Sky Reconstruction Functions
#'
#' @references \insertAllCited{}
#'
#' @return An object from the class _list_ that includes the following: (1) the
#'   reconstructed sky ([SpatRaster-class]), (2) the output produced by
#'   [fit_cie_sky_model()], (3) an object from the class _list_ that includes an
#'   object from the class `lm` (see [stats::lm()]) and the RMSE, both being the
#'   result of validating (2) with a k-fold approach and following
#'   \insertCite{Pineiro2008;textual}{rcaiman},(4) the `dist_to_plant` argument
#'   used when [fit_cie_sky_model()] was called, (5) the `sky_points` argument
#'   used when [extract_rl()] was called,  and (6) the out-of-range index (see
#'   details).
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
#' plot(bin)
#'
#'
#' mx <- optim_normalize(caim, bin)
#' caim <- normalize(caim, mx = mx, force_range = TRUE)
#' ecaim <- enhance_caim(caim, m, HSV(239, 0.85, 0.5))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#'
#' set.seed(7)
#' sky <- ootb_sky_reconstruction(r, z, a, m, bin,
#'                                sor_filter_cv = TRUE, sor_filter_dn = TRUE)
#'
#' sky$sky
#' plot(sky$sky)
#' sky$model_validation$rmse
#' plot(r/sky$sky>1.15)
#' plot(sky$model_validation$reg$model$x, sky$model_validation$reg$model$y)
#' abline(0,1)
#' error <- sky$model_validation$reg$model$x - sky$model_validation$reg$model$y
#' plot(sky$sky_points$z, error,
#'                      xlab = "zenith angle", ylab = "relative radiance error")
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#'
#' }
ootb_sky_reconstruction <- function(r, z, a, m, bin, g = NULL,
                                    sor_filter_cv = FALSE,
                                    sor_filter_dn = FALSE,
                                    refine_sun_coord = FALSE,
                                    interpolate = TRUE,
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

  .get_sky_cie <- function(model) {
    sky_cie <- cie_sky_model_raster(z, a,
                                    model$sun_coord$zenith_azimuth,
                                    model$coef) * model$zenith_dn
    names(sky_cie) <- "CIE sky"
    sky_cie
  }

  .get_metric <- function(model) {
    median(abs(model$pred - model$obs))
  }

  if (is.null(g)) {
    g <- sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
  }

  # START optimize dist_to_plant ####
  g30 <- sky_grid_segmentation(z, a, 30)
  g30[!m] <- 0
  dist_to_plant <- 11
  sampling_pct <- 0
  while (sampling_pct < 100 & dist_to_plant > 3) {
    dist_to_plant <- dist_to_plant-2
    sky_points <- extract_sky_points(r, bin, g, dist_to_plant = dist_to_plant)
    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
      (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
  }
  if (sampling_pct < 75) {
    dist_to_plant <- 1
  }
  # END optimize dist_to_plant ####

  sun_coord <- extract_sun_coord(r, z, a, bin, g)
  sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
  if (!is.null(input_sky_points)) {
    sky_points <- rbind(input_sky_points, sky_points)
  }
  rl <- extract_rl(r, z, a, sky_points, use_window = dist_to_plant != 1)
  if (sd(rl$sky_points$rl) < 0.01) {
    warning("Overexposed image")
  }

  method <- c("BFGS", "CG", "SANN")
  models <- Map(function(x) fit_cie_sky_model(rl, sun_coord, method = x),
                method)
  metric <- Map(.get_metric, models)
  i <- which.min(metric)
  model <- models[[i]]
  method <- method[i]
  sun_coord <- model$sun_coord

  # START sampling on the out-of-range zone ####
  if (dist_to_plant != 1) {
    model.2 <- model
    model$pred <- model$obs * 1e+10
    rl.2 <- rl
    while (.get_metric(model) - .get_metric(model.2) > 0.0001) {
      model <- model.2
      rl <- rl.2

      ratio <- r / .get_sky_cie(model)
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
                           use_window = dist_to_plant != 1)
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
  if (!is.null(input_sky_points)) {
    rl$sky_points <- rbind(rl$sky_points[1:nrow(input_sky_points)],
      .filter(rl$sky_points[(nrow(input_sky_points)+1):
                              nrow(rl$sky_points)], c("col", "row"), 3)
      )
  } else {
    rl$sky_points <- .filter(rl$sky_points, c("col", "row"), 3)
  }

  .sor_filter <- function(sky_points, r, fun_ct, fun_d, k, s) {
    central_tendency <- numeric(nrow(sky_points))
    dispersion <- numeric(nrow(sky_points))
    for (i in 1:nrow(sky_points)) {
      spatial_distances <- sqrt((sky_points$row - sky_points[i, "row"])^2 +
                                  (sky_points$col - sky_points[i, "col"])^2)
      phi1 <- sky_points$a * pi/180
      phi2 <- sky_points[i, "a"] * pi/180
      theta1 <- sky_points$z * pi/180
      theta2 <- sky_points[i, "z"] * pi/180
      spherical_distance <- ncol(r)/2 *
        acos(sin(phi1) * sin(phi2) + cos(phi1) *
               cos(phi2) * cos(theta1 - theta2)) %>% suppressWarnings()
      spherical_distance[is.nan(spherical_distance)] <-
        spatial_distances[is.nan(spherical_distance)]
      ds <- extract_dn(r, sky_points[, c("row", "col")])
      # u <- ds[order(spatial_distances), ][2:k + 1, 3]
      u <- ds[order(spherical_distance), ][2:k + 1, 3]
      central_tendency[i] <- fun_ct(u, na.rm = TRUE)
      dispersion[i] <- fun_d(u, na.rm = TRUE)
    }
    abs(ds[,3] - central_tendency) < (dispersion * s)
  }
  if (sor_filter_cv == TRUE) {
    cv <- terra::focal(r, 3, sd) / terra::focal(r, 3, mean)
    u <- .sor_filter(rl$sky_points, cv, mean, sd, 20, 1.5)
    if (!is.null(input_sky_points)) {
      u[1:nrow(input_sky_points)] <- TRUE
    }
    rl$sky_points <- rl$sky_points[u, ]
  }
  if (sor_filter_dn == TRUE) {
    u <- .sor_filter(rl$sky_points, r, median, IQR, 5, 1.5)
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
    x <- c(x, extract_dn(.get_sky_cie(model.2)/rl$zenith_dn,
                         rl$sky_points[folds[[i]], c("row", "col")],
                         use_window = dist_to_plant != 1)[,3])
    y <- c(y, extract_dn(r/rl$zenith_dn,
                         rl$sky_points[folds[[i]], c("row", "col")],
                         use_window = dist_to_plant != 1)[,3])
  }
  reg <- lm(x~y) #following Pineiro2008 10.1016/j.ecolmodel.2008.05.006
  # END K-fold method ####

  sky_cie <- .get_sky_cie(model)
  error <- median(abs(1 - reg$model$y/reg$model$x))
  # START interpolate ####
  if (interpolate) {
      w <-  1 - stats::plogis(error, 0.042, 0.05) + stats::plogis(0, 0.042, 0.05)
      r2 <- summary(reg) %>% .$r.squared
      if (!is.numeric(r2)) r2 <- 0
      k <- 2 + round(r2 * 10)
      p <- abs(model$coef[1]) + log(model$coef[3]+1)
      p <- 2 + sqrt(p)
      sky <- interpolate_sky_points(rl$sky_points, r, k = k, p = p,
                                    rmax = ncol(r)/7) * rl$zenith_dn
      sky <- sky * (1 - w) + sky_cie * w
      sky <- terra::cover(sky, sky_cie)
      names(sky) <- paste0("Weighted average, ",
                           "w=", round(w, 2), ", ",
                           "k=", k, ", ",
                           "p=", round(p, 2))
  } else {
    sky <- sky_cie
  }
  # END interpolate ####

  .calc_oor_index <- function(sky) {
    ratio <- r / sky
    ratio[is.infinite(ratio)] <- 1e+10
    out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
    out.of.range_ratio <- sum(out.of.range_ratio[]^2,
                              na.rm = TRUE)
    out.of.range_ratio
  }

  list(sky = sky,
       model = model,
       model_validation = list(reg = reg,
                               rmse = .calc_rmse(reg$model$y - reg$model$x)),
       dist_to_plant = dist_to_plant,
       sky_points = rl$sky_points,
       oor_index = .calc_oor_index(sky)
       )
}
