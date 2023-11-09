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
#'   the out-of-range index (see [fit_cie_sky_model()]), (4) sky points that
#'   were not involved in obtaining (2), (5) an object from the class `lm` (see
#'   [stats::lm()]) that is the result of validating (1) with (4) and the method
#'   recommended by \insertCite{Pineiro2008;textual}{rcaiman}, and (6) a
#'   binarized image produced with (1), the coefficients from (4) and
#'   [thr_mblt()] with [apply_thr()], using 'w=0.95'. If a filling source is
#'   provided, only a reconstructed sky ([SpatRaster-class]) is returned.
#'
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' m <- !is.na(z)
#' mn <- quantile(caim$Blue[m], 0.01)
#' mx <- quantile(caim$Blue[m], 0.99)
#' r <- normalize(caim$Blue, mn, mx, TRUE)
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
#' bin <- ootb_obia(caim, z, a, m, HSV(239, 0.85, 0.5), gamma = NULL)
#' plot(bin)
#' bin <- ootb_mblt(r, z, a, bin)$bin
#' plot(bin)
#'
#' set.seed(7)
#' sky <- ootb_sky_reconstruction(r, z, a, bin)
#'
#' sky$sky
#' sky$validation %>% summary()
#' plot(sky$sky)
#' plot(r/sky$sky)
#' hist(r/sky$sky, xlim = c(0, 2), breaks = 255)
#' hist((r/sky$sky)[bin], xlim = c(0, 2), breaks = 255)
#' plot((r/sky$sky)>1.1)
#'
#' plot(sky$bin)
#'
#' sky2 <- ootb_sky_reconstruction(r, z, a, sky$bin, sky$sky)
#' plot(sky2)
#' plot(r/sky2)
#' hist(r/sky2, xlim = c(0, 2), breaks = 255)
#' hist((r/sky2)[sky$bin], xlim = c(0, 2), breaks = 255)
#' plot((r/sky2)>1.1)
#' }
ootb_sky_reconstruction <- function(r, z, a, bin,
                                    filling_source = NULL,
                                    dist_to_plant = 3,
                                    sun_coord = NULL,
                                    general_sky_type = NULL,
                                    twilight = TRUE,
                                    rmse = TRUE,
                                    method = "BFGS",
                                    try_grids = TRUE,
                                    thin_points = TRUE,
                                    refine_sun_coord = TRUE,
                                    try_optims = TRUE,
                                    force_sampling = TRUE,
                                    interpolate = TRUE
                                    ) {
  if (is.null(filling_source)) {
    g <- sky_grid_segmentation(z, a, 10)
    sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
    rl <- extract_rl(r, z, a, sky_points)

    .get_r2 <- function(model) {
      if (suppressWarnings(is.null(model))) {
        r2 <- 0
      } else {
        fit <- try(lm(model$pred~model$obs), silent = TRUE)
        if (is(fit)[1] == "lm") {
          if (coefficients(fit)[2] > 0 | sd(model$pred) == 0) {
            r2 <- suppressWarnings(summary(fit)) %>% .$r.squared
          } else {
            r2 <- 0
          }
        } else {
          r2 <- 0
        }
      }
      r2
    }

    if (!is.null(sun_coord)) {
      model <- fit_cie_sky_model(r, z, a,
                                 rl$sky_points,
                                 rl$zenith_dn,
                                 rmse = rmse,
                                 sun_coord,
                                 general_sky_type = general_sky_type,
                                 twilight = twilight,
                                 method = method)
    } else {
      if (sd(rl$sky_points$rl) < 0.01) {
        sun_coord <- extract_sun_coord(r, z, a, bin,
                                       sky_grid_segmentation(z, a, 10))
        model <- fit_cie_sky_model(r, z, a,
                                   rl$sky_points,
                                   rl$zenith_dn,
                                   rmse = FALSE,
                                   sun_coord,
                                   general_sky_type = general_sky_type,
                                   twilight = twilight,
                                   method = method)
      } else {
        # ==========================
        # START sun coord extraction
        # ==========================
        m <- list(
          m_full = !is.na(z),
          m_12oclock = mask_hs(a, 0, 90) | mask_hs(a, 270, 360),
          m_3oclock = mask_hs(a, 180, 360),
          m_6oclock = mask_hs(a, 90, 270),
          m_9oclock = mask_hs(a, 0, 180)
        )
        .run_fit <- function(m, rl) {
          fit_cie_sky_model(r, z, a,
                            rl$sky_points,
                            rl$zenith_dn,
                            rmse = rmse,
                            extract_sun_coord(r, z, a, bin & m, g),
                            general_sky_type = general_sky_type,
                            twilight = twilight,
                            method = method)
        }
        models <- Map(.run_fit, m, Map(function(i) rl, 1:5))
        r2 <- Map(.get_r2, models)
        model <- models[[which.max(r2)]]
        sun_coord <- model$sun_coord
        # ========================
        # END sun coord extraction
        # ========================
      }
      # =======================
      # START Try several grids
      # =======================
      if (try_grids) {
        general_sky_type = strsplit(model$sky_type, ",")[[1]][1]
        .chessboard <- function(size) {
          cb <- chessboard(r, size)
          cb[is.na(z)] <- 0
          cb
        }
        ges <- list(
          sky_grid_segmentation(z, a, 10),
          sky_grid_segmentation(z, a, 7.5),
          sky_grid_segmentation(z, a, 5),
          .chessboard(round(ncol(r)/30)),
          .chessboard(round(ncol(r)/20))
        )
        .get_rl <- function(g, dist_to_plant) {
          sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
          extract_rl(r, z, a, sky_points)
        }

        dist_to_plant_x_2 <- dist_to_plant * 2
        if (.is_even(dist_to_plant_x_2)) {
          dist_to_plant_x_2 <- dist_to_plant_x_2 + 1
        }
        dist_to_plant_x_3 <- dist_to_plant * 3
        if (.is_even(dist_to_plant_x_3)) {
          dist_to_plant_x_3 <- dist_to_plant_x_3 + 1
        }
        dists <- list(dist_to_plant,
                      dist_to_plant_x_2,
                      dist_to_plant_x_3,
                      dist_to_plant_x_2,
                      dist_to_plant
                      )

        rls <- Map(.get_rl, ges, dists)
        models <- Map(.run_fit, m, rls)
        models[[length(models)+1]] <- model
        r2 <- Map(.get_r2, models)
        i <- which.max(r2)
        if (i != (length(models)+1)) {
          model <- models[[i]]
          rl <- rls[[i]]
          g <- ges[[i]]
          dist_to_plant <- dists[[i]]
        }
      }
      # =================================================
      # START model improvement by the thinning of points
      # =================================================
      if (thin_points) {
        general_sky_type = strsplit(model$sky_type, ",")[[1]][1]
        .improve_model <- function() {
          r2 <- .get_r2(model)
          error <- model$obs - model$pred
          i <- stats::kmeans(error, 2)
          i <- which.min(i$centers %>% abs()) == i$cluster
          # if (length(i) < no) i <- 1:length(error)
          model2 <- fit_cie_sky_model(r, z, a,
                                      rl$sky_points[i,],
                                      rmse = rmse,
                                      rl$zenith_dn,
                                      sun_coord,
                                      general_sky_type = general_sky_type,
                                      twilight = FALSE,
                                      method = method)
          if (!(is.na(model$coef) %>% any())) {
            if (r2 > .get_r2(model2)) {
              return(FALSE)
            } else {
              rl$sky_points <<- rl$sky_points[i,]
              model <<- model2
              return(TRUE)
            }
          } else {
            return(FALSE)
          }
        }

        no <- length(rl$sky_points$dn) / 2
        while(.improve_model() &
              .get_r2(model) < 0.85 &
              length(rl$sky_points$dn) > no ) NA
      }
      # ===============================================
      # START model improvement by sun coord refinement
      # ===============================================
      if (refine_sun_coord) {
        path <- system.file("external", package = "rcaiman")
        skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))
        skies <- paste0(skies$general_sky_type, ", ", skies$description)
        std_sky_no <- match(model$sky_type, skies)

        lock <- TRUE
        nn <- expand.grid(-1:1, -1:1)[-5,]
        .get_label <- function(za) paste0(za[1], "_", za[1])
        label <- "45_45"
        .get_za <- function(label) strsplit(label, "_")[[1]] %>% as.numeric()
        .zas <- c()

        while(lock) {
          row_col <- row_col_from_zenith_azimuth(z, a,
                                                 sun_coord$zenith_azimuth)$row_col
          l <- list()
          .sun_coord <- sun_coord
          za <- sun_coord$zenith_azimuth
          for (i in 1:nrow(nn)) {
            .sun_coord$zenith_azimuth[1] <- za[1] + nn[i,1]
            .sun_coord$zenith_azimuth[2] <- za[2] + nn[i,2]
            label <- .get_label(.sun_coord$zenith_azimuth)
            if (any(match(.zas, label))) {
              l[[i]] <- .sun_coord
              .zas <- c(.zas, .get_label(.sun_coord$zenith_azimuth))
            } else {
              l[[i]] <- NULL
            }
          }

          .run_fit.2 <- function(sun_coord) {
            if (suppressWarnings(is.null(sun_coord))) {
              model <- NULL
            } else {
              model <- fit_cie_sky_model(r, z, a,
                                         rl$sky_points,
                                         rl$zenith_dn,
                                         rmse = rmse,
                                         sun_coord,
                                         std_sky_no = std_sky_no,
                                         twilight = FALSE,
                                         method = method)
            }
            model
          }
          models <- Map(.run_fit.2, l)
          r2 <- Map(.get_r2, models)
          lock <- any(.get_r2(model) < r2)
          if(lock) {
            model <- models[[which.max(r2)]]
            sun_coord <- model$sun_coord
          }
        }
      }
    }
    # =============================================
    # END model improvement by sun coord refinement
    # =============================================
    .calc_oor_index <- function(sky) {
      ratio <- r / sky
      ratio[is.infinite(ratio)] <- 10000
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
      out.of.range_ratio <- sum(out.of.range_ratio[]^2,
                                          na.rm = TRUE)
      out.of.range_ratio
    }

    if (length(model$coef) != 5) {
      sky_cie <- terra::rast(z)
      names(sky_cie) <- "Fail"
    } else {
      sky_cie <- cie_sky_model_raster(z, a,
                                      model$sun_coord$zenith_azimuth,
                                      model$coef) * model$zenith_dn
      names(sky_cie) <- "CIE sky"
    }
    # =======================================================
    # START model improvement by testing optimization methods
    # =======================================================
    if (try_optims) {
      path <- system.file("external", package = "rcaiman")
      skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))
      skies <- paste0(skies$general_sky_type, ", ", skies$description)
      std_sky_no <- match(model$sky_type, skies)
      .fun <- function(i) {
        .model <- NULL
        try(.model <- fit_cie_sky_model(r, z, a,
                                        rl$sky_points,
                                        rl$zenith_dn,
                                        model$sun_coord,
                                        std_sky_no = std_sky_no,
                                        rmse = rmse,
                                        method = i))
        if (suppressWarnings(is.null(.model) | length(.model$coef) != 5)) {
          sky_cie <- terra::rast(z)
          sky_cie[] <- 1e-10
        } else {
          sky_cie <- cie_sky_model_raster(z, a,
                                          .model$sun_coord$zenith_azimuth,
                                          .model$coef) * .model$zenith_dn
        }
        list(sky_cie = sky_cie, model = .model)
      }
      l <- Map(.fun, c(rep("SANN", 10), "BFGS", "CG"))
      l[[length(l)+1]] <- list(
                sky_cie = cie_sky_model_raster(z, a,
                                               model$sun_coord$zenith_azimuth,
                                               model$coef) * model$zenith_dn,
                model = model
                )
      i <- Map(function(x) .calc_oor_index(x$sky_cie), l) %>% which.min()
      sky_cie <- l[[i]]
      model <- sky_cie$model
      sky_cie <- sky_cie$sky_cie
      names(sky_cie) <- "CIE sky"
    }
    # =============================================================
    # START model improvement by sampling on the areas out-of-range
    # =============================================================
    sky_points.2 <- NULL
    if (force_sampling) {
      ratio <- r / sky_cie
      ratio[is.infinite(ratio)] <- 10000
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
      sky_points.2 <- extract_sky_points(out.of.range_ratio, bin, g,
                                         dist_to_plant)
      sky_points.2 <- rbind(rl$sky_points[ , c("row", "col")], sky_points.2)
      rl.2 <- extract_rl(r, z, a, sky_points.2, g, no_of_points = NULL)
      rl.2$sky_points$rl <- rl.2$sky_points$rl/rl$zenith_dn
      rl.2$zenith_dn <- rl$zenith_dn

      .fun <- function(i) {
        .model <- NULL
        try(.model <- fit_cie_sky_model(r, z, a,
                                        rl.2$sky_points,
                                        rl$zenith_dn,
                                        model$sun_coord,
                                        custom_sky_coef = model$coef,
                                        rmse = rmse,
                                        method = i))

        if (suppressWarnings(is.null(.model) | length(.model$coef) != 5)) {
          sky_cie <- terra::rast(z)
          sky_cie[] <- 1e-10
        } else {
          sky_cie <- cie_sky_model_raster(z, a,
                                          .model$sun_coord$zenith_azimuth,
                                          .model$coef) * .model$zenith_dn
        }
        list(sky_cie = sky_cie, model = .model)
      }
      l <- Map(.fun, c(rep("SANN", 10), "BFGS", "CG"))
      l[[length(l)+1]] <- list(
        sky_cie = cie_sky_model_raster(z, a,
                                       model$sun_coord$zenith_azimuth,
                                       model$coef) * model$zenith_dn,
        model = model
      )
      i <- Map(function(x) .calc_oor_index(x$sky_cie), l) %>% which.min()
      sky_cie <- l[[i]]
      model <- sky_cie$model
      sky_cie <- sky_cie$sky_cie
      names(sky_cie) <- "CIE sky"
      if (i != length(l)) rl <- rl.2
    }
    # ===========================================================
    # END model improvement by sampling on the areas out-of-range
    # ===========================================================
    # =================
    # START interpolate
    # =================
    if (interpolate) {
      sky <- interpolate_sky_points(rl$sky_points, r, k = 3,
                                    rmax = ncol(r)/7)
      r2 <- .get_r2(model)
      sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
      sky_w_k3 <- terra::cover(sky, sky_cie)
      names(sky_w_k3) <- "Weighted average (k = 3)"

      sky <- interpolate_sky_points(rl$sky_points, r, k = 10,
                                    rmax = ncol(r)/7)
      sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
      sky_w_k10 <- terra::cover(sky, sky_cie)
      names(sky_w_k10) <- "Weighted average (k = 10)"

      pred <- extract_dn(sky_cie, rl$sky_points[,c("row", "col")])[,3]
      residual <-  rl$sky_points$dn - pred
      residual[abs(residual) > 10] <- 0
      residual <- cbind(rl$sky_points, residual)

      r_residual <- interpolate_sky_points(residual, r, k = 3,
                                           rmax = ncol(r)/7,
                                           col_id = "residual")
      sky <- sky_cie + r_residual
      sky_r_k3 <- terra::cover(sky, sky_cie)
      names(sky_r_k3) <- "Residual interpolation (k = 3)"

      r_residual <- interpolate_sky_points(residual, r, k = 10,
                                           rmax = ncol(r)/7,
                                           col_id = "residual")
      sky <- sky_cie + r_residual
      sky_r_k10 <- terra::cover(sky, sky_cie)
      names(sky_r_k10) <- "Residual interpolation (k = 10)"

      l <- list(sky_cie, sky_w_k3, sky_w_k10, sky_r_k3, sky_r_k10)
    } else {
      l <- list(sky_cie)
    }
    # ===============
    # END interpolate
    # ===============

    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    v <- terra::buffer(v, 2)
    bin[v] <- 0

    if (!is.null(sky_points.2)) {
      v <- cellFromRowCol(r, sky_points.2$row, sky_points.2$col) %>%
        xyFromCell(r, .) %>% vect()
      v <- terra::buffer(v, 2)
      bin[v] <- 0
    }

    bin <- bin == 1
    g <- sky_grid_segmentation(z, a, 5)
    sky_points <- extract_sky_points(r, bin, g)

    .get_rmse <- function(sky) {
      x <- extract_dn(sky, sky_points)[,3]
      y <- extract_dn(r, sky_points)[,3]
      .calc_rmse(x - y)
    }

    i <- Map(.get_rmse, l) %>% which.min()
    sky <- l[[i]]

    x <- extract_dn(sky, sky_points)[,3]
    y <- extract_dn(r, sky_points)[,3]
    fit <- lm(x~y)
    mblt <- coefficients(fit)
    r2 <- summary(fit)$r.squared
    bin <- apply_thr(r, thr_mblt(sky, mblt[1]*255, mblt[2] * r2 * 0.95))

    sky <- list(sky = sky,
                model = model,
                oor = .calc_oor_index(sky),
                sky_points = sky_points,
                validation = fit,
                bin = bin)
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
