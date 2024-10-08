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
#' @param interpolate Logical vector of length one. If `TRUE`,
#'   [interpolate_sky_points()] will be used.
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
#' sky <- ootb_sky_reconstruction(r, z, a, m, bin)
#'
#' sky$sky
#' plot(sky$sky)
#' sky$model_validation$rmse
#' plot(r/sky$sky>1.15)
#' plot(sky$model_validation$reg$model$x, sky$model_validation$reg$model$y)
#' abline(0,1)
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#'
#' # masking the sunlit canopy and running again
#' .a <- azimuth_image(z, sky$model$sun_coord$zenith_azimuth[2])
#' m_fuzzy <- (normalize((abs((.a-180)))^1.8) + normalize(sqrt(90-z))) / 2
#' sun_theta <- sky$model$sun_coord$zenith_azimuth[1]
#' m_fuzzy <- normalize(m_fuzzy, 0, (90-sun_theta)/90, TRUE)
#'
#' sky_cie <- cie_sky_model_raster(z, a, sky$model$sun_coord$zenith_azimuth,
#'                                 sky$model$coef)
#' .bin <- !apply_thr(sky_cie, thr_isodata(sky_cie[m])) & bin
#' .caim <- read_caim()
#' .mx <- optim_normalize(.caim, .bin)
#' .caim <- normalize(.caim, mx = .mx, force_range = TRUE)
#'
#' mem <- membership_to_color(.caim, HSV(239, 0.85, 0.5))
#' mem <- mem$membership_to_target_color
#' mem <- normalize(sqrt(mem*0.1), 0, 1)
#' m_fuzzy <- m * m_fuzzy + mem * (1 - m_fuzzy)
#' bin <- apply_thr(m_fuzzy*ecaim, thr_isodata(ecaim[m]))
#'
#' set.seed(7)
#' g <- sky_grid_segmentation(z, a, 5)
#' sky <- ootb_sky_reconstruction(r, z, a, m, bin, g, refine_sun_coord = TRUE)
#'
#' sky$sky
#' plot(sky$sky)
#' sky$model_validation$rmse
#' plot(r/sky$sky>1.15)
#' plot(sky$model_validation$reg$model$x, sky$model_validation$reg$model$y)
#' abline(0,1)
#'
#' plot(bin)
#' points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)
#' }
ootb_sky_reconstruction <- function(r, z, a, m, bin, g = NULL,
                                    refine_sun_coord = FALSE,
                                    interpolate = TRUE
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

  .get_custom_coef_matrix <- function(n, w) {
    custom_sky_coef <- matrix(model$coef + .noise(w), ncol = 5)
    for (i in 1:n) {
      custom_sky_coef <- rbind(custom_sky_coef,
                               model$coef + .noise(w))
    }
  }

  .get_sky_cie <- function(model) {
    sky_cie <- cie_sky_model_raster(z, a,
                                    model$sun_coord$zenith_azimuth,
                                    model$coef) * model$zenith_dn
    names(sky_cie) <- "CIE sky"
    sky_cie
  }

  .get_metric <- function(model) {
    .calc_rmse(model$pred - model$obs)
  }

  if (is.null(g)) {
    g <- sky_grid_segmentation(z, a, 10)
    i <- grep("001", g[])
    g[i] <- 1000
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

  # START sampling on the areas out-of-range ####
  if (dist_to_plant != 1) {
    # .bin <- bin
    model.2 <- model
    model$pred <- model$obs * 1e+10
    rl.2 <- rl
    while (.get_metric(model) - .get_metric(model.2) > 0.0001) {
      model <- model.2
      rl <- rl.2

      ratio <- r / .get_sky_cie(model)
      ratio[is.infinite(ratio)] <- 1e+10
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)

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
    # bin <- .bin
  }
  # END sampling on the areas out-of-range ####

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
  rmse <- .calc_rmse(1 - reg$model$x/reg$model$y)

  # START interpolate ####
  if (interpolate) {
      w <-  1 - (rmse / (0.15 + rmse)) #Michaelis Menten
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
       model_validation = list(reg = reg, rmse = rmse),
       dist_to_plant = dist_to_plant,
       sky_points = rl$sky_points,
       oor_index = .calc_oor_index(sky)
       )
}
