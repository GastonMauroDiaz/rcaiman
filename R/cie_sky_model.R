#' CIE sky model
#'
#' Written by Gaston Mauro Diaz based on Pascal code by Mait Lang.
#'
#' Angles should be provided in radians.
#'
#' @param AzP Numeric vector. Azimuth angle of a sky point.
#' @param Zp Numeric vector. Zenith Angle of a sky point.
#' @param AzS Numeric vector of length one. Azimuth angle of the sun.
#' @param Zs Numeric vector of length one. Zenith angle of the sun.
#' @param .a,.b,.c,.d,.e Numeric vector of length one. Sky model parameter.
#'
#' @noRd
#' @references http://dx.doi.org/10.1016/j.energy.2016.02.054
#'
#' @return Numeric vector of length equal to AzP length.
.cie_sky_model <- function(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e) {
  # calculate angular distance between sky point and Sun
  Chi <- .calc_angular_distance(Zp, AzP, Zs, AzS)

  # Gradation function
  Phi_Z <- 1 + .a * exp(.b / cos(Zp))
  Phi_0 <- 1 + .a * exp(.b)
  gradation <- Phi_Z / Phi_0

  # Indicatrix function
  F_Chi <-  1 + .c * (exp(.d * Chi) - exp(.d * pi/2)) + .e * cos(Chi)^2
  F_Zs  <-  1 + .c * (exp(.d *  Zs) - exp(.d * pi/2)) + .e * cos(Zs)^2
  indicatrix <- F_Chi / F_Zs

  unname(gradation * indicatrix)
}


#' CIE sky model raster
#'
#' @inheritParams sky_grid_segmentation
#' @param sun_coord Numeric vector of length 2. Zenith and azimuth angles in
#'   degrees, corresponding to the location of the solar disk center.
#' @param sky_coef Numeric vector of length 5. Parameters of the sky model.
#'
#' @family  cie sky model functions
#' @export
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1400, lens())
#' a <- azimuth_image(z)
#' path <- system.file("external", package = "rcaiman")
#' skies <- read.csv(file.path(path, "15_CIE_standard_skies.csv"))
#' # parameters are from http://dx.doi.org/10.1016/j.energy.2016.02.054
#' sky_coef <- skies[4,1:5]
#' sun_coord <- c(45, 0)
#' plot(cie_sky_model_raster(z, a, sun_coord, sky_coef))
#' }
cie_sky_model_raster <- function(z, a, sun_coord, sky_coef) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(.get_max(z) < 90)

  stopifnot(length(sun_coord) == 2)
  stopifnot(length(sky_coef) == 5)

  Zp <- degree2radian(z[])
  AzP <- degree2radian(a[])

  Zs <- degree2radian(sun_coord[1])
  AzS <- degree2radian(sun_coord[2])

  relative_luminance <- .cie_sky_model(AzP, Zp, AzS, Zs,
                                       as.numeric(sky_coef[1]),
                                       as.numeric(sky_coef[2]),
                                       as.numeric(sky_coef[3]),
                                       as.numeric(sky_coef[4]),
                                       as.numeric(sky_coef[5]))
  z[] <- relative_luminance
  z
}


#' Fit CIE sky model
#'
#' Fit CIE sky model using data from real hemispherical photographs.
#'
#' This function use maximum likelihood to estimate the parameters of the CIE
#' sky model that fit best to data sampled from a real image. It will try to
#' estimate the parameters  The result include the output produced by
#' \code{\link[bbmle]{mle2}}, the root mean square error, the adjusted r
#' squared, the digital number at the zenith, the sun coordinates (zenith and
#' azimuth angle in degrees), and the description of the standard sky from which
#' the initial parameters were drawn \insertCite{Li2016}{rcaiman}.
#'
#' @inheritParams ootb_mblt
#' @param sky_marks data.frame. The result of a call to
#'   \code{\link{extract_sky_marks}}.
#' @inheritParams cie_sky_model_raster
#' @param std_sky_no Numeric vector. Standard sky number from Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman}.
#' @param general_sky_type Character vector of length one. "Overcast", "Clear",
#'   or "Partly cloudy". See Table 1 from \insertCite{Li2016;textual}{rcaiman}.
#' @param use_window Logical vector of length one. If \code{TRUE}, a 3 by 3
#'   window will be used to extract the sky digital number from \code{r}.
#' @inheritParams bbmle::mle2
#'
#' @references \insertRef{Li2016}{rcaiman}
#'
#' @family  cie sky model functions
#'
#' @export
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' blue <- gbc(r$Blue)
#' bin <- ootb_cie_mblt(blue, z, a)$bin
#' sky_marks <- extract_sky_marks(blue, bin, g)
#' sun_coord <- extract_sun_mark(blue, bin, z, a, g)
#'
#' # since extract_sun_mark does not "see" below the horizon, this script
#' # explore til the end of civic twilight
#' sun_coords <- data.frame(z = c(sun_coord[1], seq(90, 96, 1) ),
#'                          a = sun_coord[2])
#' models <- Map(function(i) fit_cie_sky_model(blue, z, a, sky_marks,
#'                                             as.numeric(sun_coords[i,]),
#'                                             general_sky_type = "Clear"),
#'               1:nrow(sun_coords))
#' rmse <- Map(function(i) models[[i]]$rmse, seq_along(models)) %>% unlist()
#' model <- models[[which.min(rmse)]]
#' sky <- cie_sky_model_raster(z, a, model$sun_coord, model$fit@coef[-6]) *
#'   model$zenith_dn
#' plot(sky)
#' }
fit_cie_sky_model <- function(r, z, a, sky_marks, sun_coord,
                              std_sky_no = NULL,
                              general_sky_type = NULL ,
                              use_window = TRUE,
                              method = "BFGS") {
  if (!requireNamespace("bbmle", quietly = TRUE)) {
    stop(paste("Package \"bbmle\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }

  stopifnot(ncol(sky_marks) == 2)
  stopifnot(length(sun_coord) == 2)
  .check_if_r_z_and_a_are_ok(r, z, a)

  sun_a_z <- degree2radian(sun_coord[c(2,1)])
  AzS <- sun_a_z[1]
  Zs <- sun_a_z[2]

  cells <- cellFromRowCol(a, sky_marks$row, sky_marks$col)
  sky_marks$a <- a[cells]
  sky_marks$z <- z[cells]

  if (use_window) {
    xy <-  xyFromCell(r, cells)
    sky_marks$dn <-  extract(r, xy, buffer = 1.5, fun = "mean")
  } else {
    sky_marks$dn <-  r[cells]
  }

  AzP <- degree2radian(sky_marks$a)
  Zp <- degree2radian(sky_marks$z)

  z_thr <- 2
  zenith_dn <- c()
  while (length(zenith_dn) < 20) {
    zenith_dn <- sky_marks$dn[sky_marks$z < z_thr]
    z_thr <- z_thr + 2
  }
  if ((z_thr - 2) > 20) warning(paste0("Zenith DN were estimated from pure ",
                                 "sky pixels from zenith angle up to ",
                                 z_thr, "."))
  zenith_dn <- mean(zenith_dn)
  attr(zenith_dn, "max_zenith_angle") <- z_thr
  sky_marks$dn <- sky_marks$dn / zenith_dn

  flog <- function(.a, .b, .c, .d, .e, S) {
    media <- .cie_sky_model(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e)
    - sum(stats::dnorm(sky_marks$dn, mean = media, sd = exp(S)))
  }


  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  calc_rmse <- function(x) sqrt(mean(x^2))
  calc_r2 <- function(pred, obs){
    model <- lm(pred ~ obs)
    r2 <- suppressWarnings(summary(model)$adj.r.squared)
    if (is.nan(r2)) browser()
    r2
  }

  fun <- function(i) {
    fit <- NA
    try(
      fit <- bbmle::mle2(flog, list(.a = as.numeric(skies[i,1]),
                                    .b = as.numeric(skies[i,2]),
                                    .c = as.numeric(skies[i,3]),
                                    .d = as.numeric(skies[i,4]),
                                    .e = as.numeric(skies[i,5]),
                                    S = 0), method = method),
      silent = TRUE
    )

    if (any(fit@details$convergence, is.na(fit))) {
      return(NA)
    } else {
      pred <- .cie_sky_model(AzP, Zp, AzS, Zs,
                             .a = as.numeric(skies[i,1]),
                             .b = as.numeric(skies[i,2]),
                             .c = as.numeric(skies[i,3]),
                             .d = as.numeric(skies[i,4]),
                             .e = as.numeric(skies[i,5]))
      error <- pred - sky_marks$dn
      rmse <- calc_rmse(error)
      r2 <- calc_r2(pred, sky_marks$dn)
      return(list(fit = fit,
                  pred = pred,
                  obs = sky_marks$dn,
                  rmse = rmse,
                  r2 = r2,
                  zenith_dn = zenith_dn,
                  sun_coord = round(radian2degree(c(Zs, AzS))),
                  sky_type = paste0(skies[i,"general_sky_type"], ", ",
                                   skies[i,"description"])))
    }
  }

  if (!is.null(std_sky_no)) {
    skies <- skies[std_sky_no,]
  }
  if (!is.null(general_sky_type) & is.null(std_sky_no)) {
    indices <- skies$general_sky_type == general_sky_type
    skies <- skies[indices,]
  }

  fit <- suppressWarnings(Map(fun, 1:nrow(skies)))
  rmse <- Map(function(i) fit[[i]][4], 1:length(fit)) %>% unlist()
  if (all(is.na(rmse))) {
    return(list(fit = NA, pred = NA, obs = NA, rmse = NA,
                r2 = NA, zenith_dn = zenith_dn,
                sun_coord = round(radian2degree(c(Zs, AzS))),
                sky_type = NA))
  } else {
    return(fit[[which.min(rmse)]])
  }
}


#' Choose a standard CIE sky model
#'
#' Choose one from the 15th standard CIE sky models based on sky DN values.
#'
#' This function assume the same than \code{\link{find_sky_dns}}, so please
#' refer to that function help, first paragraph under Details headline.
#'
#' To kwow more about the CIE sky models, please refer to
#' \insertCite{Li2016;textual}{rcaiman}.
#'
#'
#' @inheritParams fit_cie_sky_model
#' @inheritParams extract_sky_marks
#'
#' @references \insertRef{Li2016}{rcaiman}
#'
#' @family  cie sky model functions
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' sky <- choose_std_cie_sky_model(blue, z, a)$relative_luminance
#' plot(sky)
#' }
choose_std_cie_sky_model <- function(r, z, a, bin) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  compareRaster(r, bin)

  g <- sky_grid_segmentation(z, a, 30)
  sky_marks <- extract_sky_marks(r, bin, g)
  sun_coord <- extract_sun_mark(r, bin, z, a, g)

  cells <- cellFromRowCol(a, sky_marks$row, sky_marks$col)
  #sky_marks$a <- a[cells]
  sky_marks$z <- z[cells]

  xy <-  xyFromCell(r, cells)
  sky_marks$dn <-  extract(r, xy, buffer = 1.5, fun = "mean")

  z_thr <- 2
  zenith_dn <- c()
  while (length(zenith_dn) < 20) {
    zenith_dn <- sky_marks$dn[sky_marks$z < z_thr]
    z_thr <- z_thr + 2
  }

  zenith_dn <- mean(zenith_dn)
  attr(zenith_dn, "max_zenith_angle") <- z_thr
  sky_marks$dn <- sky_marks$dn / zenith_dn

  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  calc_rmse <- function(x) sqrt(mean(x^2))
  calc_r2 <- function(pred, obs){
    model <- lm(pred ~ obs)
    suppressWarnings(summary(model)$adj.r.squared)
  }

  fun <- function(i) {
    std_sky <- cie_sky_model_raster(z, a, sun_coord, as.numeric(skies[i,1:5]))

    xy <-  xyFromCell(std_sky, cells)
    pred <-  extract(std_sky, xy, buffer = 1.5, fun = "mean")

    error <- pred - sky_marks$dn
    rmse <- calc_rmse(error)
    r2 <- calc_r2(pred, sky_marks$dn)
    return(list(relative_luminance = std_sky,
                zenith_dn = zenith_dn,
                rmse = rmse,
                r2 = r2,
                sun_coord = sun_coord,
                std_sky_no = i,
                general_sky_type = skies[i,"general_sky_type"])
           )

  }

  std_skies <- Map(fun, 1:nrow(skies))
  rmse <- Map(function(i) std_skies[[i]][2], 1:length(std_skies)) %>% unlist()
  std_skies[[which.min(rmse)]]
}



#' Autofit CIE sky model
#'
#' Autofit CIE sky model using data from real hemispherical photographs.
#'
#' By using \code{\link{ootb_cie_mblt}}, it obtains a binarized image. Then use
#' it to extract the sky marks and sun coordinates required by
#' \code{\link{fit_CIE_sky_model}} (a sky grid of 10 degrees
#' --\code{sky_grid_segmentation(z, a, 10)}-- and the default parameters of
#' \code{\link{extract_sun_marks}} and \code{\link{extract_sun_mark}} are used).
#'
#' However, the binarized image is not used as is, the \code{w} parameter (see
#' this concept on the \code{\link{thr_image}} documentation) is estimated based
#' on the sky border length calculated from the binarized image obtained with
#' \code{w} equal to \code{0.5} (by means of a Laplacian filter). The ratio of
#' the numer of pixels on the border to the total number of pixels is added to
#' \code{0.5}, but values above \code{0.9} are not allowed.
#'
#' The argument \code{general_sky_type} of \code{\link{fit_CIE_sky_model}} is
#' obtained with \code{\link{choose_std_cie_sky_model}}.
#'
#' If the general sky type is "Clear", since \code{\link{extract_sun_mark}} is
#' incapable of providing the sun position below the horizon, the sun below the
#' horizon (91 degrees) to the end of the astronomical twilight is explorer with
#' 3 degrees resolution, and the best model is returned.
#'
#' Since each general sky type has several standard CIE sky models, each of
#' those are tested and the best model is returned.
#'
#' The best model is the one with low root mean square error.
#'
#' @inheritParams fit_cie_sky_model
#'
#' @family  cie sky model functions
#'
#' @export
#'
#' @examples
#' #' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' fit <- autofit_cie_sky_model(blue, z, a)
#' plot(fit$cie_sky)
#' }
autofit_cie_sky_model <- function(r, z, a) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  .estimate_w <- function(bin) {
    sky_border <- focal(bin, matrix(c(0, 1, 0,
                                      1,-4, 1,
                                      0, 1, 0), nrow = 3))
    no_of_px_on_circle <- sum(!is.na(bin)[])
    no_of_px_on_sky_border <- sum((sky_border != 0)[], na.rm = TRUE)
    w <- 0.5 + no_of_px_on_sky_border / no_of_px_on_circle
    if (w > 0.9) w <- 0.9
    w
  }

  .autofit <- function(r, z, a,
                       sky_marks,
                       sun_coord,
                       general_sky_type) {
    if (general_sky_type == "Clear") {
      sun_coords <- data.frame(z = c(sun_coord[1], seq(90, 112, 2)),
                               a = sun_coord[2])
      models <- Map(function(i) {
                      suppressWarnings(
                        fit_cie_sky_model(r, z, a, sky_marks,
                                          as.numeric(sun_coords[i, ]),
                                          general_sky_type = general_sky_type,
                                          use_window = TRUE)
                      )
                    },
                    1:nrow(sun_coords))
      rmse <- Map(function(i) models[[i]]$rmse, seq_along(models)) %>% unlist()
      model <- models[[which.min(rmse)]]
    } else {
      model <- suppressWarnings(
                fit_cie_sky_model(r, z, a, sky_marks,
                                  sun_coord,
                                  general_sky_type = general_sky_type,
                                  use_window = TRUE)
               )
    }
    model
  }

  mblt <- ootb_cie_mblt(r, z, a)
  w <- .estimate_w(mblt$bin)
  bin <- r / mblt$sky > w
  g <- sky_grid_segmentation(z, a, 10)
  sky_marks <- extract_sky_marks(r, bin, g)
  sun_coord <- extract_sun_mark(r, mblt$bin, z, a, g)

  # l <- list()
  # l$overcast <- .autofit(r, z, a, sky_marks, sun_coord, "Overcast")
  # l$partly <- .autofit(r, z, a, sky_marks, sun_coord, "Partly cloudy")
  # l$clear <- .autofit(r, z, a, sky_marks, sun_coord, "Clear")
  #
  # .fun <- function(x) sum(x[x > 1 | x < 0]^2)
  # error <- c(.fun(l$overcast$lr), .fun(l$partly$lr), .fun(l$clear$lr))
  # model <- l[[which.min(error)]]

  model <- .autofit(r, z, a, sky_marks, sun_coord, mblt$general_sky_type)

  if (!is.na(model$r2)) {
    model$relative_luminance <- cie_sky_model_raster(z, a,
                                           model$sun_coord,
                                           model$fit@coef[-6]) * model$zenith_dn
  }

  model$mblt <- mblt
  model$w <- w
  model$sky_marks <- sky_marks

  no_cells <- index <- round(z) == model$sun_coord[1] &
                                                round(a) == model$sun_coord[2]
  no_cells[] <- 1:ncell(index)
  sun_mark <- rowColFromCell(index, no_cells[index][1])
  model$sun_mark <- sun_mark
  model
}

