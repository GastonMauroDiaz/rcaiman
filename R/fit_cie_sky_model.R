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
#' @inheritParams ootb_mblt
#' @param sun_coord Numeric vector of length two. Zenith and azimuth angles in
#'   degrees, corresponding to the location of the solar disk center.
#' @param sky_coef Numeric vector of length five. Parameters of the sky model.
#'
#' @family  Sky Reconstruction Functions
#'
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
  .is_single_layer_raster(z)
  .is_single_layer_raster(a)
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  stopifnot(length(sun_coord) == 2)
  stopifnot(length(sky_coef) == 5)

  Zp <- .degree2radian(z[])
  AzP <- .degree2radian(a[])

  Zs <- .degree2radian(sun_coord[1])
  AzS <- .degree2radian(sun_coord[2])

  relative_luminance <- .cie_sky_model(AzP, Zp, AzS, Zs,
                                       as.numeric(sky_coef[1]),
                                       as.numeric(sky_coef[2]),
                                       as.numeric(sky_coef[3]),
                                       as.numeric(sky_coef[4]),
                                       as.numeric(sky_coef[5]))
  terra::values(z) <- relative_luminance
  z
}



#' Fit CIE sky model
#'
#' Use maximum likelihood to estimate the coefficients of the CIE sky model that
#' best fit to data sampled from a real scene.
#'
#' This function assumes an hemispherical image as input. It is based on
#' \insertCite{Lang2010;textual}{rcaiman}. In theory, the best result would be
#' obtained with data showing a linear relation between digital numbers and the
#' amount of light reaching the sensor. However, because the CIE sky model is
#' indeed the adjoin of two mathematical model, it is capable of handling any
#' non-linearity since it is not a physical model with strict assumptions.
#'
#' Ultimately, when the goal is to calculate the ratio of canopy to sky digital
#' numbers, if the latter is accurately constructed, any non-linearity will be
#' canceled. Please, see \code{\link{interpolate_sky_points}} for further
#' considerations.
#'
#' Nevertheless, the recommended input for this function is data pre-processed
#' with the HSP software package \insertCite{Lang2013}{rcaiman}. Please, refer
#' to \code{\link{write_sky_points}} for additional details about HSP.
#'
#' The following code exemplifies how this package can be used to compare the
#' manually-guided fitting provided by HSP against the automatic fitting
#' provided by this package. The code assumes that the user is working within an
#' RStudio project located in the HSP project folder.
#'
#' \preformatted{
#' r <- read_caim("manipulate/IMG_1013.pgm") %>% normalize()
#' plot(r)
#' z <- zenith_image(ncol(r), lens())
#' a <- azimuth_image(z)
#' manual_input <- read_manual_input(".", "IMG_1013" )
#' sun_coord <- manual_input$sun_mark$row_col
#' sun_coord <- zenith_azimuth_from_row_col(r, sun_coord, lens())
#' sky_points <- manual_input$sky_marks
#' rl <- extract_rl(r, z, a, sky_points)
#' model <- fit_cie_sky_model(r, z, a, rl$sky_points, rl$zenith_dn, sun_coord)
#' cie_sky <- model$relative_luminance * model$zenith_dn
#' plot(r/cie_sky)
#'
#' sky_coef <- read_opt_sky_coef(".", "IMG_1013")
#' cie_sky_manual <- cie_sky_model_raster(z, a, sun_coord$zenith_azimuth, sky_coef)
#' img <- read_caim("manipulate/IMG_1013.pgm")
#' cie_sky_manual <- cie_sky_manual * manual_input$zenith_dn
#' plot(img/cie_sky_manual)
#' }
#'
#' If you use this function in your research, please cite
#' \insertCite{Lang2010;textual}{rcaiman} in addition to this package.
#'
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_coneshaped_model
#' @param zenith_dn Numeric vector of length 1. Zenith digital number, see
#'   \code{\link{extract_rl}} for how to obtain it.
#' @param sun_coord An object of class \emph{list}. The result of a call to
#'   \code{\link{extract_sun_coord}}.
#' @inheritParams cie_sky_model_raster
#' @param std_sky_no Numeric vector. Standard sky number from Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman}.
#' @param general_sky_type Character vector of length one. It could be any of
#'   these: "Overcast", "Clear", or "Partly cloudy". See Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman} for additional details.
#' @param twilight Logical vector of length one. If it is \code{TRUE} and the
#'   initial standard parameters belong to the "Clear" general sky type, sun
#'   zenith angles from 90 to 96 degrees will be tested (civic twilight). This
#'   is necessary since \code{\link{extract_sun_coord}} would mistakenly
#'   recognize the center of what can be seen of the solar corona as the solar
#'   disk.
#' @param rmse Logical vector of length one. If it is \code{TRUE}, the criteria
#'   for selecting the best sky model is to choose the one with less root mean
#'   square error calculated by using \code{sky_points} as reference values.
#'   Otherwise, the criteria is to evaluate the whole hemisphere by calculating
#'   the product between the square ratio of \code{r} to the sky model and the
#'   fraction of pixels from this new layer that are above one or below zero,
#'   and selecting the sky model that produce the least value.
#' @inheritParams bbmle::mle2
#'
#' @references \insertAllCited{}
#'
#' @return The result includes the following: (1) the output produced by
#'   \code{\link[bbmle]{mle2}}, (2) the 5 coefficients, (3) observed and
#'   predicted values, the sun coordinates --zenith and azimuth angle in
#'   degrees--, (4) the relative luminance image calculated for every pixel
#'   using the estimated coefficients and corresponding sun coordinates, (4) the
#'   digital number at the zenith, and (5) the description of the standard sky
#'   from which the initial coefficients were drawn. See
#'   \insertCite{Li2016;textual}{rcaiman} to know more about these coefficients.
#'
#' @family  Sky Reconstruction Functions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' sky_blue_sample <- crop(caim, ext(610,643,760,806))
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' ecaim <- enhance_caim(caim, !is.na(z), sky_blue, gamma = 2.2)
#' bin <- apply_thr(ecaim, 0.75)
#' g <- sky_grid_segmentation(z, a, 10)
#' r <- gbc(caim$Blue*255)
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' sky_points <- extract_sky_points(r, bin, g)
#' rl <- extract_rl(r, z, a, sky_points)
#' model <- fit_cie_sky_model(r, z, a, rl$sky_points,
#'                            rl$zenith_dn, sun_coord,
#'                            rmse = TRUE,
#'                            general_sky_type = "Partly cloudy")
#' sky_cie <- model$relative_luminance * model$zenith_dn
#' sky_cie <- normalize(sky_cie, 0, 1, TRUE)
#' plot(sky_cie)
#' plot(r/sky_cie)
#' }
fit_cie_sky_model <- function(r, z, a, sky_points, zenith_dn, sun_coord,
                              std_sky_no = NULL,
                              general_sky_type = NULL ,
                              twilight = TRUE,
                              rmse = FALSE,
                              method = "BFGS") {
  if (!requireNamespace("bbmle", quietly = TRUE)) {
    stop(paste("Package \"bbmle\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }

  stopifnot(is.data.frame(sky_points))
  stopifnot(length(sun_coord$zenith_azimuth) == 2)
  .check_if_r_z_and_a_are_ok(r, z, a)

  AzP <- .degree2radian(sky_points$a)
  Zp <- .degree2radian(sky_points$z)

  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  .fun <- function(i) {
    sun_a_z <- .degree2radian(rev(sun_coord$zenith_azimuth))
    AzS <- sun_a_z[1]
    Zs <- sun_a_z[2]

    flog <- function(.a, .b, .c, .d, .e, S) {
      media <- .cie_sky_model(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e)
      - sum(stats::dnorm(sky_points$rl, mean = media, sd = exp(S)))
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
                  sun_coord = sun_coord,
                  sky_type = NA))
    } else {
      pred <- .cie_sky_model(AzP, Zp, AzS, Zs,
                             .a = as.numeric(skies[i,1]),
                             .b = as.numeric(skies[i,2]),
                             .c = as.numeric(skies[i,3]),
                             .d = as.numeric(skies[i,4]),
                             .e = as.numeric(skies[i,5]))
      relative_luminance <- cie_sky_model_raster(z, a,
                                                 sun_coord$zenith_azimuth,
                                                 fit@coef[-6])
      return(list(mle2_output = fit,
                  coef = fit@coef[-6],
                  obs = sky_points$rl,
                  pred = pred,
                  relative_luminance = relative_luminance,
                  zenith_dn = zenith_dn,
                  sun_coord = sun_coord,
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
        sun_coord$zenith_azimuth <-  c(civic_twilight[i],
                                      sun_coord$zenith_azimuth[2])
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
      sun_rl <- sky[x$sun_coord$row_col[1], x$sun_coord$row_col[2]][,]
      if (x$sun_coord$zenith_azimuth[1] < 90 &
          sun_rl > quantile(sky[], 0.9, na.rm = TRUE)
          ) {
        ratio <- r / sky
        ratio[is.infinite(ratio)] <- -0.1
        m <- ratio < 0 | ratio > 1
        area_outside_expected_values <- sum(m[], na.rm = TRUE)
        w <- area_outside_expected_values / total_area
        return(sum(ratio[]^2, na.rm = TRUE) * w)
      } else {
        if (x$sun_coord$zenith_azimuth[1] >= 90) {
          ratio <- r / sky
          ratio[is.infinite(ratio)] <- -0.1
          m <- ratio < 0 | ratio > 1
          area_outside_expected_values <- sum(m[], na.rm = TRUE)
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
    if (model$sun_coord$zenith_azimuth[1] >= 90) {
        warning(paste("Sun zenith angle was overwriten so \"row_col\"",
                      "and \"zenith_azimuth\" does not match as in the",
                      "original \"sun_coord\" argument.",
                      "If you need to recalculate the \"sun_coord\",",
                      "check the row_col_from_zenith_azimuth() function."))
      }
    model$relative_luminance[is.infinite(model$relative_luminance)] <- 0
    return(model)
  }
}

