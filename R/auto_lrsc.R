auto_lrsc <- function(r, z, a, general_sky_type = NULL) {
  .auto_lrsc <- function(r, z, a, general_sky_type = "Clear") {
    g <- sky_grid_segmentation(z, a, 10)

    mblt <- ootb_mblt(r, z, a)

    sky_marks <- extract_sky_marks(r, mblt$bin, g,
                                   dist_to_plant = 3,
                                   min_raster_dist = 6
    ) #
    sun_coord <- extract_sun_mark(r, mblt$bin, z, a, g)

    if (general_sky_type == "Clear") {
      sun_coords <- data.frame(
        z = c(sun_coord[1], seq(91, 118, 3)), # 118 astro-dusk
        a = sun_coord[2]
      )
      models <- Map(
        function(i) {
          fit_cie_sky_model(r, z, a, sky_marks,
                            as.numeric(sun_coords[i, ]),
                            general_sky_type = general_sky_type,
                            use_kernel = TRUE
          )
        },
        1:nrow(sun_coords)
      )
      rmse <- Map(function(i) models[[i]]$rmse, seq_along(models)) %>% unlist()
      model <- models[[which.min(rmse)]]
    } else {
      model <- fit_cie_sky_model(r, z, a, sky_marks,
                                 sun_coord,
                                 general_sky_type = general_sky_type,
                                 use_kernel = TRUE
      )
    }

    cie_sky <- cie_sky_model_raster(z, a, model$sun_coord, model$fit@coef[-6]) *
      model$zenith_dn

    w <- model$r2
    sky <- cie_sky * w + mblt$sky * (1 - w)

    lr <- r / sky
    list(lr = lr, sky = sky, cie_sky = model, mblt = mblt)
  }

  if (is.null(general_sky_type)) {
    lr <- list()
    lr$overcast <- .auto_lrsc(r, z, a, "Overcast")
    lr$partly <- .auto_lrsc(r, z, a, "Partly cloudy")
    lr$clear <- .auto_lrsc(r, z, a, "Clear")

    thr <- 1.1
    error <- c(freq(lr$overcast$lr > thr, value = 1),
               freq(lr$partly$lr > thr, value = 1),
               freq(lr$clear$lr > thr, value = 1))
    lr <- lr[[which.min(error)]]

  } else {
    lr <- .auto_lrsc(r, z, a, general_sky_type)
  }
  lr
}
