.auto_lrsc <- function(r, z, a,
                       mblt,
                       sky_marks,
                       sun_coord,
                       general_sky_type,
                       use_window = TRUE) {
  if (general_sky_type == "Clear") {
    sun_coords <- data.frame(z = c(sun_coord[1], seq(91, 118, 3)), # 118 astro
                             a = sun_coord[2])
    models <- Map(function(i) {
      fit_cie_sky_model(r, z, a, sky_marks,
                        as.numeric(sun_coords[i, ]),
                        general_sky_type = general_sky_type,
                        use_window = use_window)
    },
    1:nrow(sun_coords))
    rmse <- Map(function(i) models[[i]]$rmse, seq_along(models)) %>% unlist()
    model <- models[[which.min(rmse)]]
  } else {
    model <- fit_cie_sky_model(r, z, a, sky_marks,
                               sun_coord,
                               general_sky_type = general_sky_type,
                               use_window = use_window)
  }

  if (is.na(model$r2)) {
    sky <- mblt$sky
  } else {
    cie_sky <- cie_sky_model_raster(z, a,
                                    model$sun_coord,
                                    model$fit@coef[-6]) * model$zenith_dn
    w <- model$r2
    sky <- cie_sky * w + mblt$sky * (1 - w)
  }

  lr <- r / sky
  list(lr = lr, sky = sky,
       cie_sky = model,
       mblt = mblt)
}


auto_lrsc <- function(r, z, a, general_sky_type = NULL) {
  .check_if_r_z_and_a_are_ok(r, z, a)


  mblt <- cie_mblt(r, z, a)

  # estimate w
  sky_border <- focal(mblt$bin, matrix(c(0, 1, 0,
                                         1,-4, 1,
                                         0, 1, 0), nrow = 3))
  no_of_px_on_circle <- sum(!is.na(z)[])
  no_of_px_on_sky_border <- sum((sky_border != 0)[], na.rm = TRUE)
  w <- 0.5 + no_of_px_on_sky_border / no_of_px_on_circle
  if (w > 0.9) w <- 0.9

  bin <- r / mblt$sky > w

  g <- sky_grid_segmentation(z, a, 10)
  sky_marks <- extract_sky_marks(r, bin, g,
                                 dist_to_plant = 3,
                                 min_raster_dist = 3)
  sun_coord <- extract_sun_mark(r, mblt$bin, z, a, g)

  if (is.null(general_sky_type)) {
    l <- list()
    l$overcast <- .auto_lrsc(r, z, a, mblt, sky_marks, sun_coord, "Overcast")
    l$partly <- .auto_lrsc(r, z, a, mblt, sky_marks, sun_coord, "Partly cloudy")
    l$clear <- .auto_lrsc(r, z, a, mblt, sky_marks, sun_coord, "Clear")

    .fun <- function(x) sum(x[x > 1 | x < 0]^2)
    error <- c(.fun(l$overcast$lr), .fun(l$partly$lr), .fun(l$clear$lr))
    lr <- l[[which.min(error)]]

  } else {
    lr <- .auto_lrsc(r, z, a, mblt, general_sky_type)
  }
  lr$g <- g
  lr$bin <- bin
  lr
}
