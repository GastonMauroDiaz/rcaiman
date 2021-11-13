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
#' This function use maximum likelihood to estimate the coefficients of the CIE
#' sky model that fit best to data sampled from a real image. The result include
#' the output produced by \code{\link[bbmle]{mle2}}, the 5 model coefficients,
#' observed and predicted values, the sun coordinates (zenith and azimuth angle
#' in degrees), the relative luminance calculated for every pixel using the
#' estimated coefficients and corresponding sun coordinate, the digital number
#' at the zenith, and the description of the standard sky from which the initial
#' coefficients were drawn. See \insertCite{Li2016;textual}{rcaiman} to known
#' more about these coefficients.
#'
#' This function assume an hemispherical image as input. It is based on
#' \insertCite{Lang2010;textual}{rcaiman}. In theory, the better result would be
#' obtained with data showing a linear relation between digital numbers and the
#' amount of light reaching the sensor. However, the CIE sky model is indeed the
#' adjoin of two mathematical models, one controlling the gradation between the
#' zenith and the horizon (two parameters), and the other controlling the
#' gradation originated at the solar disk (three parameters). This make the CIE
#' model capable of cope with any non-linearity.
#'
#' Ultimately, if the goal is to calculate the ratio of canopy DN to sky DN, if
#' the latter is accurately constructed, any non-linearity will be canceled.
#' Please, see \code{\link{interpolate_dns}} for further considerations.
#'
#' @inheritParams ootb_mblt
#' @param sky_marks An object of class data.frame. The result of a call to
#'   \code{\link{extract_sky_marks}}.
#' @param sun_mark An object of class list. The result of a call to
#'   \code{\link{extract_sun_mark}}.
#' @inheritParams cie_sky_model_raster
#' @param std_sky_no Numeric vector. Standard sky number from Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman}.
#' @param general_sky_type Character vector of length one. "Overcast", "Clear",
#'   or "Partly cloudy". See Table 1 from \insertCite{Li2016;textual}{rcaiman}.
#' @param use_window Logical vector of length one. If \code{TRUE}, a 3 by 3
#'   window will be used to extract the sky digital number from \code{r}.
#' @param twilight Logical vector of length one. If it is \code{TRUE} and the
#'   initial standard parameters belong to the "Clear" general sky type, sun
#'   zenith angles from 90 to 96 degrees will be tested (civic twilight). This
#'   is necessary since \code{\link{extract_sun_mark}} would mistakenly
#'   recognize the gravity center of what can be seen of the solar corona as the
#'   solar disk.
#' @param rmse Logical vector of length one. If it is \code{TRUE}, the criteria
#'   for selecting the best sky model is to choose the one with lest root mean
#'   square error calculated from the sample (sky marks). Otherwise, the
#'   criteria is to evaluate the whole hemisphere by calculating the ratio of
#'   \code{r} to the sky model, then square it, and selecting the sky model that
#'   produce the least value.
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
#' bin <- ootb_mblt(blue, z, a, is_horizon_visible = TRUE)$bin
#' sky_marks <- extract_sky_marks(blue, bin, g)
#' sun_mark <- extract_sun_mark(blue, bin, z, a, g)
#' model <- fit_cie_sky_model(blue, z, a, sky_marks, sun_mark)
#' plot(model$relative_luminance)
#' }
fit_cie_sky_model <- function(r, z, a, sky_marks, sun_mark,
                              std_sky_no = NULL,
                              general_sky_type = NULL ,
                              use_window = TRUE,
                              twilight = TRUE,
                              rmse = FALSE,
                              method = "BFGS") {
  if (!requireNamespace("bbmle", quietly = TRUE)) {
    stop(paste("Package \"bbmle\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }

  stopifnot(ncol(sky_marks) == 2)
  stopifnot(length(sun_mark$zenith_azimuth) == 2)
  .check_if_r_z_and_a_are_ok(r, z, a)

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

  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  .fun <- function(i) {
    sun_a_z <- degree2radian(rev(sun_mark$zenith_azimuth))
    AzS <- sun_a_z[1]
    Zs <- sun_a_z[2]

    flog <- function(.a, .b, .c, .d, .e, S) {
      media <- .cie_sky_model(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e)
      - sum(stats::dnorm(sky_marks$dn, mean = media, sd = exp(S)))
    }

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

    if (any(try(fit@details$convergence, silent = TRUE), is.na(fit))) {
      return(list(mle2_output = NA,
                  coef = NA,
                  obs = NA,
                  pred = NA,
                  relative_luminance = NA,
                  zenith_dn = zenith_dn,
                  sun_mark = sun_mark,
                  sky_type = NA))
    } else {
      pred <- .cie_sky_model(AzP, Zp, AzS, Zs,
                             .a = as.numeric(skies[i,1]),
                             .b = as.numeric(skies[i,2]),
                             .c = as.numeric(skies[i,3]),
                             .d = as.numeric(skies[i,4]),
                             .e = as.numeric(skies[i,5]))
      relative_luminance <- cie_sky_model_raster(z, a,
                                                 sun_mark$zenith_azimuth,
                                                 fit@coef[-6])
      return(list(mle2_output = fit,
                  coef = fit@coef[-6],
                  obs = sky_marks$dn,
                  pred = pred,
                  relative_luminance = relative_luminance,
                  zenith_dn = zenith_dn,
                  sun_mark = sun_mark,
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

  fit <- suppressWarnings(Map(.fun, 1:nrow(skies)))

  if (twilight) {
    indices <- match(11:15, as.numeric(rownames(skies)))
    indices <- indices[!is.na(indices)]
    if (length(indices) != 0) {
      skies <- skies[indices,]
      civic_twilight <-  c(seq(90, 96, 1))
      for (i in seq_along(civic_twilight)) {
        sun_mark$zenith_azimuth <-  c(civic_twilight[i],
                                      sun_mark$zenith_azimuth[2])
        fit <- c(fit, suppressWarnings(Map(.fun, 1:nrow(skies))))
      }
    }
  }

  total_area <- sum(!is.na(z)[], na.rm = TRUE)
  .calc_ratio_squared <- function(x) {
    if (length(x$coef) != 5) {
      return(NA)
    } else {
      sky <- x$relative_luminance * x$zenith_dn
      sun <- sky[x$sun_mark$row_col[1], x$sun_mark$row_col[2]]
      if (x$sun_mark$zenith_azimuth[1] < 90 &
          sun > quantile(sky[], 0.9, na.rm = TRUE)
          ) {
        m <- sky < 0 | sky > 1
        area_outside_expected_values <- sum(m[], na.rm = TRUE)
        ratio <- r / sky
        w <- area_outside_expected_values / total_area
        return(sum(ratio[]^2, na.rm = TRUE) * w)
      } else {
        if (x$sun_mark$zenith_azimuth[1] >= 90) {
          m <- sky < 0 | sky > 1
          area_outside_expected_values <- sum(m[], na.rm = TRUE)
          ratio <- r / sky
          w <- area_outside_expected_values / total_area
          return(sum(ratio[]^2, na.rm = TRUE) * w)
        } else {
          return(NA)
        }
      }
    }
  }

  if (!rmse) {
    error <- Map(.calc_ratio_squared, fit) %>% unlist()
  } else {
    error <- Map(function(x) .calc_rmse(x$pred - x$obs), fit) %>% unlist()
  }

  if (all(is.na(error))) {
    return(fit[[1]])
  } else {
    model <- fit[[which.min(error)]]
    if (model$sun_mark$zenith_azimuth[1] >= 90) {
        warning(paste("Sun zenith angle was overwriten so \"row_col\"",
                      "and \"zenith_azimuth\" does not match as in the",
                      "original \"sun_mark\" argument.",
                      "If you need to recalculate the \"sun_mark\",",
                      "check the row_col_from_zenith_azimuth() function."))
      }
    return(model)
  }
}

