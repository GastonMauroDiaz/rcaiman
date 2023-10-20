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
#' sky digital numbers will be available for almost every pixel, either from `r`
#' gaps or from the filling source---an exception is a filling source with
#' plenty of `NA` values, which should not be provided.
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_trend_surface
#' @inheritParams fit_cie_sky_model
#' @param dist_to_plant Numeric vector of length one or `NULL`. See
#'   [extract_sky_points()].
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
#'   binarized image produced with (1) and the coefficients from (4). If a
#'   filling source is provided, only a reconstructed sky ([SpatRaster-class])
#'   is returned.
#'
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' r <- normalize(caim$Blue)
#'
#' mn_mx <- optim_normalize(caim, !is.na(z))
#' caim <- normalize(caim, mn_mx[1], mn_mx[2], TRUE)
#'
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' sky_blue <- extract_dn(caim, sky_points, fun = median)
#'
#' bin <- ootb_obia(caim, z, a, sky_blue = sky_blue, gamma = NULL)
#' plot(bin)
#'
#' .fun <- function(method) {
#'   ootb_sky_reconstruction(r, z, a, bin, method = method)
#' }
#' set.seed(7)
#' l <- future.apply::future_Map(.fun,
#'                               c("BFGS", "CG", "SANN"), future.seed = TRUE)
#'
#' r2 <- Map(function(sky) sky$validation %>% summary() %>% .$r.squared, l) %>%
#'   unlist()
#' oor <- Map(function(sky) sky$OOR, l) %>% unlist()
#' n <- Map(function(sky) sky$sky_points %>% nrow(), l) %>% unlist()
#' i <- which.min(oor/r2/n)
#' sky <- l[[i]]
#' sky$validation %>% summary()
#'
#' sky$sky
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
                                    general_sky_type = NULL,
                                    twilight = TRUE,
                                    method = "BFGS"
                                    ) {
  if (is.null(filling_source)) {
    g <- sky_grid_segmentation(z, a, 10)
    sky_points <- extract_sky_points(r, bin, g, dist_to_plant)
    rl <- extract_rl(r, z, a, sky_points)

    if (sd(rl$sky_points$rl) < 0.01) {
      sun_coord <- extract_sun_coord(r, z, a, bin,
                                     sky_grid_segmentation(z, a, 10))
      model <- fit_cie_sky_model(r, z, a,
                                 rl$sky_points,
                                 rl$zenith_dn,
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
      .run_fit <- function(m) {
        fit_cie_sky_model(r, z, a,
                          rl$sky_points,
                          rl$zenith_dn,
                          rmse = TRUE,
                          extract_sun_coord(r, z, a, bin & m, g),
                          general_sky_type = general_sky_type,
                          twilight = twilight,
                          method = method)
      }
      models <- Map(.run_fit, m)
      .get_r2 <- function(model) {
        # plot(model$obs, model$pred)
        fit <- lm(model$pred~model$obs)
        # abline(fit)
        if (coefficients(fit)[2] > 0 | sd(model$pred) == 0) {
          r2 <- summary(fit) %>% .$r.squared
        } else {
          r2 <- 0
        }
        r2
      }
      r2 <- Map(.get_r2, models)
      model <- models[[which.max(r2)]]
      sun_coord <- model$sun_coord
      # =================================================
      # START model improvement by the thinning of points
      # =================================================
      general_sky_type = strsplit(model$sky_type, ",")[[1]][1]
      .improve_model <- function() {
        r2 <- .get_r2(model)
        error <- model$obs - model$pred
        i <- stats::kmeans(error, 2)
        i <- which.min(i$centers %>% abs()) == i$cluster
        model2 <- fit_cie_sky_model(r, z, a,
                                    rl$sky_points[i,],
                                    rmse = TRUE,
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
        }
      }

      no <- length(rl$sky_points$dn) / 2
      while(.improve_model() &
            .get_r2(model) < 0.85 &
            length(rl$sky_points$dn) > no ) NA

      # ===============================================
      # START model improvement by sun coord refinement
      # ===============================================
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
                                       rmse = TRUE,
                                       sun_coord,
                                       std_sky_no = std_sky_no,
                                       twilight = FALSE,
                                       method = method)
          }
          model
        }
        models <- Map(.run_fit.2, l)
        .get_r2 <- function(model) {
          if (suppressWarnings(is.null(model))) {
            r2 <- 0
          } else {
            # plot(model$obs, model$pred)
            fit <- lm(model$pred~model$obs)
            # abline(fit)
            if (coefficients(fit)[2] > 0 | sd(model$pred) == 0) {
              r2 <- suppressWarnings(summary(fit)) %>% .$r.squared
            } else {
              r2 <- 0
            }
          }
          r2
        }
        r2 <- Map(.get_r2, models)
        lock <- any(.get_r2(model) < r2)
        if(lock) {
          model <- models[[which.max(r2)]]
          sun_coord <- model$sun_coord
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

    # =======================================================
    # START model improvement by testing optimization methods
    # =======================================================

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
                                      rmse = TRUE,
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

    # =============================================================
    # START model improvement by sampling on the areas out-of-range
    # =============================================================
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
                                      rmse = TRUE,
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

    # ===========================================================
    # END model improvement by sampling on the areas out-of-range
    # ===========================================================
    sky <- interpolate_sky_points(rl.2$sky_points, r, k = 3,
                                  rmax = ncol(r)/7)
    r2 <- .get_r2(model)
    sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
    sky_w_k3 <- terra::cover(sky, sky_cie)
    names(sky_w_k3) <- "Weighted average (k = 3)"

    sky <- interpolate_sky_points(rl.2$sky_points, r, k = 10,
                                  rmax = ncol(r)/7)
    sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
    sky_w_k10 <- terra::cover(sky, sky_cie)
    names(sky_w_k10) <- "Weighted average (k = 10)"

    pred <- extract_dn(sky_cie, rl.2$sky_points[,c("row", "col")])[,3]
    residual <-  rl.2$sky_points$dn - pred
    residual <- cbind(rl.2$sky_points, residual)
    residual <- interpolate_sky_points(residual, r, k = 3,
                                       rmax = ncol(r)/7, col_id = "residual")
    sky <- sky_cie + residual
    sky_r_k3 <- terra::cover(sky, sky_cie)
    names(sky_r_k3) <- "Residual interpolation (k = 3)"

    residual <-  rl.2$sky_points$dn - pred
    residual <- cbind(rl.2$sky_points, residual)
    residual <- interpolate_sky_points(residual, r, k = 10,
                                       rmax = ncol(r)/7, col_id = "residual")
    sky <- sky_cie + residual
    sky_r_k10 <- terra::cover(sky, sky_cie)
    names(sky_r_k10) <- "Residual interpolation (k = 10)"

    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    v <- terra::buffer(v, 5)
    bin[v] <- 0

    v <- cellFromRowCol(r, sky_points.2$row, sky_points.2$col) %>%
      xyFromCell(r, .) %>% vect()
    v <- terra::buffer(v, 5)
    bin[v] <- 0

    bin <- bin == 1
    g <- sky_grid_segmentation(z, a, 5)
    sky_points <- extract_sky_points(r, bin, g)

    .get_rmse <- function(sky) {
      x <- extract_dn(sky, sky_points)[,3]
      y <- extract_dn(r, sky_points)[,3]
      .calc_rmse(x - y)
    }

    l <- list(sky_cie, sky_w_k3, sky_w_k10, sky_r_k3, sky_r_k10)
    i <- Map(.get_rmse, l) %>% which.min()
    sky <- l[[i]]
    # sky <- normalize(sky, 0, 1, TRUE)


    x <- extract_dn(sky, sky_points)[,3]
    y <- extract_dn(r, sky_points)[,3]
    fit <- lm(x~y)
    mblt <- coefficients(fit)
    r2 <- summary(fit)$r.squared
    bin <- apply_thr(r, thr_mblt(sky, mblt[1]*255, mblt[2] * r2 * 0.95))

    sky <- list(sky = sky,
                model = model,
                OOR = .calc_oor_index(sky),
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
