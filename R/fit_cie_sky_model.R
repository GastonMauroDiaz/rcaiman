#' Fit CIE sky model
#'
#' Use maximum likelihood to estimate the coefficients of the CIE sky model that
#' best fit to data sampled from a canopy photograph.
#'
#' This function is based on \insertCite{Lang2010;textual}{rcaiman}. In theory,
#' the best result would be obtained with data showing a linear relation between
#' digital numbers and the amount of light reaching the sensor. See
#' [extract_radiometry()] and [read_caim_raw()] for further details. As a
#' compromise solution, [gbc()] can be used.
#'
#' The following code exemplifies how this package can be used to compare the
#' manually-guided fitting provided by HSP \insertCite{Lang2013}{rcaiman}
#' against the automatic fitting provided by this package. The code assumes that
#' the user is working within an RStudio project located in the HSP project
#' folder.
#'
#' ````
#' r <- read_caim("manipulate/IMG_1013.pgm") %>% normalize()
#' z <- zenith_image(ncol(r), lens())
#' a <- azimuth_image(z)
#' manual_input <- read_manual_input(".", "IMG_1013" )
#' sun_coord <- manual_input$sun_coord$row_col
#' sun_coord <- zenith_azimuth_from_row_col(z, sun_coord, lens())
#' sky_points <- manual_input$sky_points
#' rl <- extract_rl(r, z, a, sky_points)
#' model <- fit_cie_sky_model(r, z, a, rl$sky_points, rl$zenith_dn, sun_coord)
#' cie_sky <- model$relative_luminance * model$zenith_dn
#' plot(r/cie_sky)
#'
#' r <- read_caim("manipulate/IMG_1013.pgm")
#' sky_coef <- read_opt_sky_coef(".", "IMG_1013")
#' cie_sky_manual <- cie_sky_model_raster(z, a, sun_coord$zenith_azimuth, sky_coef)
#' cie_sky_manual <- cie_sky_manual * manual_input$zenith_dn
#' plot(r/cie_sky_manual)
#' ````
#'
#'
#' @note
#'
#' If you use this function in your research, please cite
#' \insertCite{Lang2010;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`).
#'
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_coneshaped_model
#' @param zenith_dn Numeric vector of length 1. Zenith digital number, see
#'   [extract_rl()] for how to obtain it.
#' @param sun_coord An object of class *list*. The result of a call to
#'   [extract_sun_coord()] or an object with same structure and names. See also
#'   [row_col_from_zenith_azimuth()] in case you want to provide values based on
#'   date and time of acquisition and the R package 'suncalc'.
#' @inheritParams cie_sky_model_raster
#' @param custom_sky_coef Numeric vector of length five. Custom starting
#'   coefficients of the sky model. By default, they are drawn from standard
#'   skies.
#' @param std_sky_no Numeric vector. Standard sky number from Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman}.
#' @param general_sky_type Character vector of length one. It could be any of
#'   these: "Overcast", "Clear", or "Partly cloudy". See Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman} for additional details.
#' @param twilight Logical vector of length one. If it is `TRUE` and the initial
#'   standard parameters belong to the "Clear" general sky type, sun zenith
#'   angles from 90 to 96 degrees will be tested (civic twilight). This is
#'   necessary since [extract_sun_coord()] would mistakenly recognize the center
#'   of what can be seen of the solar corona as the solar disk.
#' @param rmse Logical vector of length one. If it is `TRUE`, the criteria for
#'   selecting the best sky model is to choose the one with the lowest __root
#'   mean square error (RMSE)__ calculated by using the `sky_points` argument as
#'   the source of reference values. Otherwise, the criteria is to evaluate the
#'   whole image by calculating the __out-of-range index__ as \eqn{\sum_{i =
#'   1}^{N}(r_i/sky_i)^2}, where \eqn{r} is the `r` argument, \eqn{sky} is the
#'   raster obtained from the fitted model with [cie_sky_model_raster()] and
#'   `zenith_dn`, \eqn{i} is the index that represents the position of a given
#'   pixel on the raster grid, \eqn{N} is the total number of pixels that
#'   satisfy either of these inequalities: \eqn{r_i/sky_i<0} and
#'   \eqn{r_i/sky_i>1}.
#' @inheritParams bbmle::mle2
#'
#' @references \insertAllCited{}
#'
#' @return object from the class _list_. The result includes the following: (1)
#'   the output produced by [bbmle::mle2()], (2) the 5 coefficients, (3 and 4)
#'   observed and predicted values, (5) the digital number at the zenith, (6)
#'   the sun coordinates --zenith and azimuth angle in degrees--, and (7) the
#'   description of the standard sky from which the initial coefficients were
#'   drawn. See \insertCite{Li2016;textual}{rcaiman} to know more about these
#'   coefficients.
#'
#' @family  Sky Reconstruction Functions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # Manual method after Lang et al. (2010)
#' # ImageJ can be used to digitize points
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' xy <- c(210, 451) #'originally captured with click() after x11()
#' sun_coord <- zenith_azimuth_from_row_col(z, z, c(nrow(z) - xy[2],xy[1]))
#' points(sun_coord$row_col[2], nrow(caim) - sun_coord$row_col[1],
#'        col = 3, pch = 1)
#'
#' rl <- extract_rl(caim$Blue, z, a, sky_points)
#'
#' set.seed(7)
#' model <- fit_cie_sky_model(caim$Blue, z, a, rl$sky_points,
#'                            rl$zenith_dn, sun_coord,
#'                            general_sky_type = "Clear",
#'                            rmse = FALSE,
#'                            twilight = FALSE,
#'                            method = "SANN")
#' summary(model$mle2_output)
#' plot(model$obs, model$pred)
#' abline(0,1)
#' r2 <- lm(model$pred~model$obs) %>% summary(.) %>% .$r.squared
#' r2
#' sky_cie <- cie_sky_model_raster(z, a,
#'                                 model$sun_coord$zenith_azimuth,
#'                                 model$coef) * model$zenith_dn
#' sky_cie <- normalize(sky_cie, 0, 1, TRUE)
#' plot(sky_cie)
#' plot(caim$Blue/sky_cie)
#'
#' # a quick demonstration of how to use interpolation to improve sky modelling
#' # after Lang et al. (2010)
#' sky <- interpolate_sky_points(rl$sky_points, caim$Blue, rmax = ncol(caim)/7)
#' plot(sky)
#' sky <- sky * rl$zenith_dn * (1 - r2) + sky_cie * r2
#' sky <- terra::cover(sky, sky_cie)
#' plot(sky)
#' plot(caim$Blue/sky)
#'
#' # how to provide a custom starting coefficient
#' model <- fit_cie_sky_model(caim$Blue, z, a, rl$sky_points,
#'                            rl$zenith_dn, sun_coord,
#'                            custom_sky_coef = model$coef,
#'                            method = "SANN")
#' plot(model$obs, model$pred, ylim = range(model$obs))
#' abline(0,1)
#' }
fit_cie_sky_model <- function(r, z, a, sky_points, zenith_dn, sun_coord,
                              custom_sky_coef = NULL,
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
  if (!is.null(custom_sky_coef)) {
    stopifnot(length(custom_sky_coef) == 5)
    stopifnot(is.numeric(custom_sky_coef))
    skies[1, 1:5] <- custom_sky_coef
    skies <- skies[1,]
    skies$general_sky_type <- "custom"
    skies$description <- "custom"
  }

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
      return(list(mle2_output = fit,
                  coef = fit@coef[-6],
                  obs = sky_points$rl,
                  pred = pred,
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

  # total_area <- sum(!is.na(z)[], na.rm = TRUE)
  # .calc_sasr <- function(x) {
  #   if (length(x$coef) != 5) {
  #     return(NA)
  #   } else {
  #     relative_luminance <- cie_sky_model_raster(z, a,
  #                                                x$sun_coord$zenith_azimuth,
  #                                                x$coef)
  #     sky <- relative_luminance * x$zenith_dn
  #     sun_rl <- sky[x$sun_coord$row_col[1], x$sun_coord$row_col[2]][,]
  #     if (x$sun_coord$zenith_azimuth[1] < 90 &
  #         sun_rl > quantile(sky[], 0.9, na.rm = TRUE)
  #         ) {
  #       ratio <- r / sky
  #       ratio[is.infinite(ratio)] <- -0.1
  #       m <- ratio < 0 | ratio > 1
  #       area_outside_expected_values <- sum(m[], na.rm = TRUE)
  #       w <- area_outside_expected_values / total_area
  #       return(sum(ratio[]^2, na.rm = TRUE) * w)
  #     } else {
  #       if (x$sun_coord$zenith_azimuth[1] >= 90) {
  #         ratio <- r / sky
  #         ratio[is.infinite(ratio)] <- -0.1
  #         m <- ratio < 0 | ratio > 1
  #         area_outside_expected_values <- sum(m[], na.rm = TRUE)
  #         w <- area_outside_expected_values / total_area
  #         return(sum(ratio[]^2, na.rm = TRUE) * w)
  #       } else {
  #         return(NA)
  #       }
  #     }
  #   }
  # }

  if (!rmse) {
    .calc_oor_index <- function(x) {
      if (length(x$coef) != 5) {
        sky <- terra::rast(z)
        sky[] <- 2
      } else {
        relative_luminance <- cie_sky_model_raster(z, a,
                                                   x$sun_coord$zenith_azimuth,
                                                   x$coef)
        sky <- relative_luminance * x$zenith_dn
      }

      ratio <- r / sky
      ratio[is.infinite(ratio)] <- 10000
      out.of.range_ratio <- ratio - normalize(ratio, 0, 1, TRUE)
      out.of.range_ratio <- sum(out.of.range_ratio[]^2,
                                na.rm = TRUE)
      out.of.range_ratio
    }
    error <- Map(.calc_oor_index, fit) %>% unlist()
  } else {
    error <- Map(function(x) .calc_rmse(x$pred - x$obs), fit) %>% unlist()
  }

  if (all(is.na(error))) {
    return(fit[[1]])
  } else {
    model <- fit[[which.min(error)]]
    if (model$sun_coord$zenith_azimuth[1] >= 90) {
        model$sun_coord$row_col <- c(NA, NA)
        # warning(paste("Sun zenith angle was overwriten so \"row_col\"",
        #               "and \"zenith_azimuth\" does not match as in the",
        #               "original \"sun_coord\" argument.",
        #               "If you need to recalculate the \"sun_coord\",",
        #               "check the row_col_from_zenith_azimuth() function."))
      }
    return(model)
  }
}

