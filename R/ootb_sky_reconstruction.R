#' Out-of-the-box sky reconstruction
#'
#' Build an above canopy image from a single below canopy image
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions [fit_cie_sky_model()] and [interpolate_sky_points()].
#'
#' The pipeline is an automatic version of the
#' \insertCite{Lang2010;textual}{rcaiman} method.
#'
#' Providing a `filling source` triggers an alternative pipeline in which the
#' sky is fully reconstructed with [interpolate_sky_points()] after a dense
#' sampling (\eqn{1 \times 1} degree cells), which is supported by the fact that
#' sky digital numbers will be available for every pixel, either from `r` gaps
#' or from the filling source.
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_trend_surface
#' @inheritParams fit_cie_sky_model
#' @param dist_to_plant Numeric vector of length one or `NULL`. See
#'   [extract_sky_points()].
#' @param try_grids Logical vector of length one.
#' @param thin_points Logical vector of length one.
#' @param refine_sun_coord Logical vector of length one.
#' @param try_optims Logical vector of length one.
#' @param force_sampling Logical vector of length one.
#' @param interpolate Logical vector of length one. If `TRUE`,
#'   [interpolate_sky_points()] will be used.
#'
#' @export
#'
#' @family Sky Reconstruction Functions
#'
#' @references \insertAllCited{}
#'
#' @return If a filling source is not provided, the result is an object from the
#'   class _list_ that includes the following: (1) the reconstructed sky
#'   ([SpatRaster-class]), (2) the output produced by [fit_cie_sky_model()], (3)
#'   the general sky type resulting from finding to which types belong the
#'   standard sky with the coefficients closer to the fitted coefficients (it
#'   may be different to the `general_sky_type' argument when it is provided)
#'   (4) the sky points used to fit (2), and (5) the out-of-range index (see
#'   [fit_cie_sky_model()]). If a filling source is provided, only a
#'   reconstructed sky ([SpatRaster-class]) is returned.
#'
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r <- normalize(caim$Blue)
#'
#' bin <- find_sky_pixels(r, z, a)
#' bin <- ootb_mblt(r, z, a, bin)
#' plot(bin$bin)
#'
#' mx <- optim_normalize(caim, m)
#'
#' r <- normalize(caim$Blue)
#' caim <- normalize(caim, mx = mx, force_range = TRUE)
#'
#' ecaim <- enhance_caim(caim, m, HSV(239, 0.85, 0.5))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#'
#' dist_to_plant <- optimize_dist_to_plant(r, z, a, !is.na(z), bin)
#' set.seed(7)
#' sky <- ootb_sky_reconstruction(r, z, a, bin,
#'                                dist_to_plant = dist_to_plant,
#'                                try_grids = FALSE,
#'                                thin_points = FALSE,
#'                                refine_sun_coord = FALSE,
#'                                force_sampling = FALSE,
#'                                try_optims = FALSE,
#'                                interpolate = TRUE)
#' sky$oor
#' plot(sky$sky)
#' plot(r/sky$sky)
#' plot(r/sky$sky>1.05)
#' sky$general_sky_type
#' }
ootb_sky_reconstruction <- function(r, z, a, bin,
                                    filling_source = NULL,
                                    dist_to_plant = 3,
                                    sun_coord = NULL,
                                    custom_sky_coef = NULL,
                                    std_sky_no = NULL,
                                    general_sky_type = NULL,
                                    twilight = TRUE,
                                    rmse = TRUE,
                                    method = "BFGS",
                                    try_grids = TRUE,
                                    thin_points = TRUE,
                                    refine_sun_coord = TRUE,
                                    force_sampling = TRUE,
                                    try_optims = TRUE,
                                    interpolate = TRUE
                                    ) {
  if (is.null(filling_source)) {

    if (try_grids == TRUE) {
      if (!requireNamespace("imager", quietly = TRUE)) {
        stop(paste("Package \"imager\" needed for this function to work.",
                   "Please install it."),
             call. = FALSE)
      }
    }

    if (dist_to_plant != 1) {
      size <- dist_to_plant * 2
    } else {
      size <- dist_to_plant
    }
    .this_requires_EBImage()
    kern <- EBImage::makeBrush(dist_to_plant, "box")
    dist_to_plant_img <- bin
    dist_to_plant_img <- EBImage::erode(as.array(dist_to_plant_img), kern) %>%
      terra::setValues(dist_to_plant_img, .)
    dist_to_plant_img[is.na(dist_to_plant_img)] <- 0
    dist_to_plant_img <- as.logical(dist_to_plant_img)

    g <- sky_grid_segmentation(z, a, 10)
    if (is.null(sun_coord)) sun_coord <- extract_sun_coord(r, z, a, bin, g)
    sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
    rl <- extract_rl(r, z, a, sky_points, use_window = dist_to_plant != 1)
    if (sd(rl$sky_points$rl) < 0.01) {
      warning("Overexposed image")
      try_grids = FALSE
      thin_points = FALSE
      refine_sun_coord = FALSE
      try_optims = FALSE
    }

    .get_r2 <- function(model) {
      r2 <- NA
      fit <- try(lm(model$pred~model$obs), silent = TRUE)
      if (is(fit)[1] == "lm") {
        if (coefficients(fit)[2] > 0 | sd(model$pred) == 0) {
          r2 <- suppressWarnings(summary(fit)) %>% .$r.squared
        }
      }
      if (is.na(r2)) r2 <- 0
      r2
    }

    .get_sky_cie <- function(model) {
      if (suppressWarnings(is.null(model) | length(model$coef) != 5)) {
        sky_cie <- terra::rast(z)
        sky_cie[] <- 1e-10
        names(sky_cie) <- "failed CIE sky"
      } else {
        sky_cie <- cie_sky_model_raster(z, a,
                                        model$sun_coord$zenith_azimuth,
                                        model$coef) * model$zenith_dn
        names(sky_cie) <- "CIE sky"
      }
      sky_cie
    }

    .get_cie_fit_index <- function(model, in_gaps = TRUE) {
      r2 <- .get_r2(model)
      if (r2 == 0) {
        return(1e+10)
      } else {
        if (in_gaps) {
          ratio <- r/.get_sky_cie(model)
          fit <- lm(model$pred~model$obs)
          return(1 - r2 +
                   abs(0 - coef(fit)[1]) +
                   abs(1 - coef(fit)[2]) +
                   abs(1 - median(ratio[dist_to_plant_img]))
                 )
        } else {
          fit <- lm(model$pred~model$obs)
          return(1 - r2 +
                   abs(0 - coef(fit)[1]) +
                   abs(1 - coef(fit)[2])
          )
        }
      }
    }

    .calc_oor_index <- function(sky) {
      ratio <- r / sky
      ratio[is.infinite(ratio)] <- 1e+10
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
      out.of.range_ratio <- sum(out.of.range_ratio[]^2,
                                na.rm = TRUE)
      out.of.range_ratio
    }

    path <- system.file("external", package = "rcaiman")
    skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))
    .find_std_sky_no <- function(model) {
      coef <- model$coef
      delta <- skies[, 1:5]
      for (i in 1:nrow(skies)) {
        delta[i,] <- delta[i,] - coef
      }
      apply(delta^2, 1, sum) %>% which.min()
    }

    .find_general_sky_type <- function(model) {
      std_sky_no <- .find_std_sky_no(model)
      skies[std_sky_no, "general_sky_type"]
    }

    .calc_noise <- function(general_sky_type = NULL) {
      if (is.null(general_sky_type)) {
        .sd <- apply((skies[, 1:5]), 2, sd)
      } else {
        u <- skies[,"general_sky_type"] == general_sky_type
        .sd <- apply((skies[u, 1:5]), 2, sd)
      }
      Map(function(i) rnorm(1, 0, .sd[i]), 1:5) %>% unlist()
    }

    model <- fit_cie_sky_model(r, z, a,
                               rl$sky_points,
                               rl$zenith_dn,
                               sun_coord,
                               custom_sky_coef = custom_sky_coef,
                               std_sky_no = custom_sky_coef,
                               general_sky_type = general_sky_type,
                               twilight = twilight,
                               rmse = rmse,
                               method = method)

    if (any(try_grids, thin_points, refine_sun_coord,
            force_sampling, try_optims)) {
      .model <- model
      .g <- g
      .rl <- rl
    }

    # START Try several grids ####
    if (try_grids) {
      general_sky_type <- .find_general_sky_type(model)
      .chessboard <- function(size) {
        cb <- chessboard(r, size)
        cb[is.na(z)] <- 0
        cb
      }
      ges <- list(
        sky_grid_segmentation(z, azimuth_image(z, 5), 10),
        sky_grid_segmentation(z, a, 30),
        sky_grid_segmentation(z, a, 15),
        sky_grid_segmentation(z, a, 7.5),
        .chessboard(round(ncol(r)/30)),
        .chessboard(round(ncol(r)/20)),
        .chessboard(round(ncol(r)/10))
      )
      .get_rl <- function(g) {
        sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
        extract_rl(r, z, a, sky_points, use_window = dist_to_plant != 1)
      }
      rls <- Map(.get_rl, ges)
      .run_fit <- function(m, rl) {
        fit_cie_sky_model(r, z, a,
                          rl$sky_points,
                          rl$zenith_dn,
                          sun_coord,
                          general_sky_type = general_sky_type,
                          twilight = FALSE,
                          rmse = TRUE,
                          method = method)
      }
      models <- Map(.run_fit, !is.na(z), rls)
      models[[length(models)+1]] <- model
      rls[[length(rls)+1]] <- rl
      ges[[length(ges)+1]] <- g
      cie_fit_index <- Map(.get_cie_fit_index, models)
      i <- which.min(cie_fit_index)
      model <- models[[i]]
      rl <- rls[[i]]
      g <- ges[[i]]
    }
    # END Try several grids ####

    # START model improvement by the thinning of points ####
    if (thin_points) {
      general_sky_type <- .find_general_sky_type(model)
      .improve_model <- function() {
        cie_fit_index <- .get_cie_fit_index(model)
        error <- model$obs - model$pred
        i <- stats::kmeans(error, 2)
        i <- which.min(i$centers %>% abs()) == i$cluster
        model.2 <- fit_cie_sky_model(r, z, a,
                                     rl$sky_points[i,],
                                     rl$zenith_dn,
                                     sun_coord,
                                     general_sky_type = general_sky_type,
                                     twilight = FALSE,
                                     rmse = TRUE,
                                     method = method)
        if (!(is.na(model$coef) %>% any())) {
          if (cie_fit_index < .get_cie_fit_index(model2)) {
            return(FALSE)
          } else {
            rl$sky_points <<- rl$sky_points[i,]
            model <<-model.2
            return(TRUE)
          }
        } else {
          return(FALSE)
        }
      }

      while(.improve_model() &
            length(rl$sky_points$dn) > 20 ) NA
    }
    # END model improvement by the thinning of points ####

    # START model improvement by sun_coord refinement ####
    if (refine_sun_coord) {
      std_sky_no <- .find_std_sky_no(model)
      .refine_sun_coord <- function(zenith, azimuth) {
        sun_coord <- row_col_from_zenith_azimuth(z, a, c(zenith, azimuth))
        model <- fit_cie_sky_model(r, z, a,
                                   rl$sky_points,
                                   rl$zenith_dn,
                                   sun_coord,
                                   std_sky_no = std_sky_no,
                                   twilight = FALSE,
                                   method = method)
        .get_cie_fit_index(model)
      }
      fit <- bbmle::mle2(.refine_sun_coord,
                         list(zenith = sun_coord$zenith_azimuth[1],
                              azimuth = sun_coord$zenith_azimuth[2]),
                         method = "BFGS")
      sun_coord <- row_col_from_zenith_azimuth(z, a, fit@coef)
    }
    # END model improvement by sun coord refinement ####

    sky_cie <- .get_sky_cie(model)

    # START model improvement by sampling on the areas out-of-range ####
    sky_points.2 <- NULL
    if (force_sampling) {
      ratio <- r / sky_cie
      ratio[is.infinite(ratio)] <- 1e+10
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
      sky_points.2 <- extract_sky_points(out.of.range_ratio, bin, g,
                                         dist_to_plant)
      sky_points.2 <- rbind(rl$sky_points[ , c("row", "col")], sky_points.2)
      rl.2 <- extract_rl(r, z, a, sky_points.2, g, no_of_points = NULL,
                         use_window = dist_to_plant != 1)
      rl.2$sky_points$rl <- rl.2$sky_points$rl/rl$zenith_dn
      rl.2$zenith_dn <- rl$zenith_dn

      model.2 <- NULL
      try(model2 <- fit_cie_sky_model(r, z, a,
                                      rl.2$sky_points,
                                      rl$zenith_dn,
                                      model$sun_coord,
                                      custom_sky_coef = model$coef,
                                      twilight = FALSE,
                                      method = method))

      models <- list(model,model.2)
      cie_fit_index <- Map(.get_cie_fit_index, models)
      i <- which.min(cie_fit_index)
      if (i == 2) {
        rl <- rl.2
        model <- models[[i]]
        sky_cie <- .get_sky_cie(model)
      }
    }
    # END model improvement by sampling on the areas out-of-range ####

    # START model improvement by testing optimization methods ####
    if (try_optims) {
      # general_sky_type <- skies[std_sky_no,"general_sky_type"]
      # std_sky_no <- .find_std_sky_no(model)
      .fun <- function(i) {
       model.2 <- NULL
        try(model2 <- fit_cie_sky_model(r, z, a,
                                        rl$sky_points,
                                        rl$zenith_dn,
                                        model$sun_coord,
                                        # std_sky_no = std_sky_no,
                                        custom_sky_coef = model$coef +
                                          .calc_noise(),
                                        method = i))
       model.2
      }
      .methods <- c(rep("SANN", 5), rep("BFGS",5), rep("CG",5))
      models <- Map(.fun, .methods)
      models[[length(models)+1]] <- model
      cie_fit_index <- Map(.get_cie_fit_index, models)
      i <- which.min(cie_fit_index)
      # oor <- Map(function(x) .calc_oor_index(.get_sky_cie(x)), models)
      # i <- which.min(oor)
      model <- models[[i]]
      sky_cie <- .get_sky_cie(model)
      method <- .methods[i]
    }
    # END model improvement by testing optimization methods ####

    if (any(try_grids, thin_points, refine_sun_coord,
            force_sampling, try_optims)) {
      if (.calc_oor_index(.get_sky_cie(.model)) <
          .calc_oor_index(.get_sky_cie(model))) {
        model <- .model
        #g <- .g
        rl <- .rl
      }
    }

    general_sky_type <- .find_general_sky_type(model)

    # START R^2 computation with the K-fold method ####
    folds <- seq_along(rl$sky_points$row)
    folds <- split(folds, 1:5) %>% suppressWarnings()
    r2 <- c()
    for (i in 1:5) {
      model.2 <- fit_cie_sky_model(r, z, a, rl$sky_points[-folds[[i]],],
                             rl$zenith_dn,
                             sun_coord,
                             custom_sky_coef = model$coef +
                               .calc_noise(general_sky_type),
                             twilight = FALSE,
                             method = method)
      x <- extract_dn(sky_cie, rl$sky_points[folds[[i]],1:2])[,3]
      y <- extract_dn(r, rl$sky_points[folds[[i]],1:2])[,3]
      fit <- lm(x~y)
      r2 <- c(r2,
              as.numeric(try(summary(fit)$r.squared, silent = TRUE)) %>%
                suppressWarnings())
    }
    r2 <- median(r2, na.rm = TRUE)
    # END R^2 computation with the K-fold method ####

    # START interpolate ####
    if (interpolate) {
      if (general_sky_type == "Clear" | general_sky_type == "Overcast") {
        sky <- interpolate_sky_points(rl$sky_points, r, k = 10,
                                      rmax = ncol(r)/7)
        sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
        sky <- terra::cover(sky, sky_cie)
        names(sky) <- "Weighted average (k = 10)"
      }
      if (general_sky_type == "Partly cloudy") {
        sky <- interpolate_sky_points(rl$sky_points, r, k = 3,
                                      rmax = ncol(r)/7)
        sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
        sky <- terra::cover(sky, sky_cie)
        names(sky) <- "Weighted average (k = 3)"
      }
    } else {
      sky <- sky_cie
    }
    # END interpolate ####

    sky <- list(sky = sky,
                model = model,
                general_sky_type = general_sky_type,
                sky_points = rl$sky_points,
                model_validation = fit,
                oor = .calc_oor_index(sky))
  } else {
    terra::compareGeom(r, filling_source)
    r[!bin] <- filling_source[!bin]
    sky_points <- extract_sky_points(r, !is.na(z),
                                     sky_grid_segmentation(z, a, 1),
                                     dist_to_plant = NULL,
                                     min_raster_dist = NULL)
    sky_points <- extract_rl(r, z, a, sky_points, NULL,
                             use_window = FALSE)$sky_points
    sky_points <- sky_points[!is.na(sky_points$rl),]
    sky <- interpolate_sky_points(sky_points, r)
    sky[is.na(z)] <- 0
    names(sky) <- "Reconstructed sky"
  }
  sky
}
