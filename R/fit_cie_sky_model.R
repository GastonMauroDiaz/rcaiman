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
#' model <- fit_cie_sky_model(rl, sun_coord)
#' cie_sky <- model$relative_luminance * model$zenith_dn
#' plot(r/cie_sky)
#'
#' r <- read_caim("manipulate/IMG_1013.pgm")
#' sky_coef <- read_opt_sky_coef(".", "IMG_1013")
#' cie_sky_m <- cie_sky_image(z, a, sun_coord$zenith_azimuth, sky_coef)
#' cie_sky_m <- cie_sky_manual * manual_input$zenith_dn
#' plot(r/cie_sky_m)
#' ````
#'
#'
#' @note
#'
#' The [point selection tool of ‘ImageJ’
#' software](https://imagej.net/ij/docs/guide/146-19.html#sec:Multi-point-Tool)
#' can be used to manually digitize points and create a CSV file from which to
#' read coordinates (see Examples). After digitizing the points on the image,
#' use the dropdown menu Analyze>Measure to open the Results window. To obtain
#' the CSV file, use File>Save As...
#'
#' The [QGIS software](https://qgis.org/) can also be used to manually digitize
#' points. In order to do that, drag and drop the image in an empty project,
#' create an new vector layer, digitize points manually, save the editions, and
#' close the project. To create the new vector layer go to the dropdown menu
#' Layer>Create Layer>New Geopackage Layer...
#'
#' Choose "point" in the Geometry type dropdown list and make sure the CRS is
#' EPSG:7589. To be able to input the points, remember to click first on the
#' Toogle Editing icon, and then on the Add Points Feature icon.
#'
#' If you use this function in your research, please cite
#' \insertCite{Lang2010;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`)).
#'
#'
#' @inheritParams ootb_mblt
#' @inheritParams fit_coneshaped_model
#'
#' @param rl An object of class *list*. The output of [extract_rl()] or an
#'   object with same structure and names.
#' @param sun_coord An object of class *list*. The output of
#'   [extract_sun_coord()] or an object with same structure and names. See also
#'   [row_col_from_zenith_azimuth()] in case you want to provide values based on
#'   date and time of acquisition and the `suncalc` package.
#' @inheritParams cie_sky_image
#' @param custom_sky_coef Numeric vector of length five or a numeric matrix with
#'   five columns. Custom starting coefficients of the sky model. By default,
#'   they are drawn from standard skies.
#' @param std_sky_no Numeric vector. Standard sky number from
#'   \insertCite{Li2016;textual}{rcaiman}'s Table 1.
#' @param general_sky_type Character vector of length one. It could be any of
#'   these: "Overcast", "Clear", or "Partly cloudy". See Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman} for additional details.
#' @param twilight Numeric vector of length one. Sun zenith angle (in degrees).
#'   If the sun zenith angle provided through the `sun_coord` argument is below
#'   this value, sun zenith angles from 90 to 96 degrees (civic twilight) will
#'   be tested when selecting initial optimization parameters from the general
#'   sky types _Clear_ or _Partially Cloudy_ (specifically, for standard sky
#'   numbers 7 to 15). This adjustment is necessary because
#'   [extract_sun_coord()] can mistakenly identify the visible center of the
#'   solar corona as the solar disk. Since [extract_sun_coord()] cannot output a
#'   zenith angle below 90 degrees, setting this value to 90 is equivalent to
#'   disabling this step.
#' @inheritParams bbmle::mle2
#'
#' @references \insertAllCited{}
#'
#' @return An object of the class *list*. The result includes the following: (1)
#'   the output produced by [bbmle::mle2()], (2) the 5 coefficients of the CIE
#'   model, (3) observed values, (4) predicted values, (5) the digital number at
#'   the zenith, (6) the sun coordinates (zenith and azimuth angle in degrees),
#'   (7) the optimization method (see [bbmle::mle2()]), and the initial values
#'   for optimizer (see [bbmle::mle2()]). To lear more about these initial
#'   values, see \insertCite{Li2016;textual}{rcaiman}. If [bbmle::mle2()] does
#'   not converge, (1) will be `NA` and (2) will contain the coefficients of a
#'   standard sky (the one with less RMSE when more than one is tried).
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
#' # Manual method following Lang et al. (2010)
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
#' # Idem for QGIS
#' path <- system.file("external/sky_points.gpkg",
#'                     package = "rcaiman")
#' sky_points <- terra::vect(path)
#' sky_points <- terra::extract(caim, sky_points, cells = TRUE)
#' sky_points <- terra::rowColFromCell(caim, sky_points$cell) %>% as.data.frame()
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' xy <- c(210, 451) #taken with click() after x11(), then hardcoded here
#' sun_coord <- zenith_azimuth_from_row_col(z, a, c(nrow(z) - xy[2],xy[1]))
#' points(sun_coord$row_col[2], nrow(caim) - sun_coord$row_col[1],
#'        col = 3, pch = 1)
#'
#' rl <- extract_rl(caim$Blue, z, a, sky_points)
#'
#' set.seed(7)
#' model <- fit_cie_sky_model(rl, sun_coord,
#'                            general_sky_type = "Clear",
#'                            twilight = 90,
#'                            method = "CG")
#' summary(model$mle2_output)
#' plot(model$obs, model$pred)
#' abline(0,1)
#' lm(model$pred~model$obs) %>% summary()
#'
#' sky_cie <- cie_sky_image(z, a,
#'                                 model$sun_coord$zenith_azimuth,
#'                                 model$coef) * model$zenith_dn
#' plot(sky_cie)
#' plot(caim$Blue/sky_cie)
#' }
fit_cie_sky_model <- function(rl, sun_coord,
                              custom_sky_coef = NULL,
                              std_sky_no = NULL,
                              general_sky_type = NULL ,
                              twilight = 60,
                              method = "BFGS") {
  if (!requireNamespace("bbmle", quietly = TRUE)) {
    stop(paste("Package \"bbmle\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }

  stopifnot(general_sky_type == "Overcast" |
            general_sky_type == "Partly cloudy" |
            general_sky_type == "Clear" |
            is.null(general_sky_type))
  stopifnot(is.data.frame(rl$sky_points))
  stopifnot(length(sun_coord$zenith_azimuth) == 2)

  AzP <- .degree2radian(rl$sky_points$a)
  Zp <- .degree2radian(rl$sky_points$z)

  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))
  if (!is.null(custom_sky_coef)) {
    if (is.vector(custom_sky_coef)) {
      stopifnot(length(custom_sky_coef) == 5)
      stopifnot(is.numeric(custom_sky_coef))
      skies <- skies[1,]
      skies[1, 1:5] <- custom_sky_coef
    } else {
      stopifnot(ncol(custom_sky_coef) == 5)
      custom_sky_coef <- as.matrix(custom_sky_coef)
      stopifnot(is.numeric(custom_sky_coef))
      skies <- skies[1:nrow(custom_sky_coef),]
      skies[,1:5] <- custom_sky_coef
    }
    skies$general_sky_type <- "custom"
    skies$description <- "custom"
  }

  .fun <- function(i) {
    sun_a_z <- .degree2radian(rev(sun_coord$zenith_azimuth))
    AzS <- sun_a_z[1]
    Zs <- sun_a_z[2]

    # **Note**: when the .c or .e parameters are negative, the sun can be darker
    # than the sky, which does not make sense. So, the code avoids that
    # possibility by using abs()
    # **New note**: Positive values of b create unrealistic values at the
    # horizon. Positive values of d produce negative values, which are
    # unrealistc
    flog <- function(.a, .b, .c, .d, .e, S) {
      x <- .cie_sky_model(AzP, Zp, AzS, Zs,
                          .a, -abs(.b), abs(.c), -abs(.d), abs(.e))
      stats::median(abs(x - rl$sky_points$rl)) #because some points might be wrong
      # .calc_rmse(x - rl$sky_points$rl)
      # .calc_rmse(1 - rl$sky_points$rl/x )
    }

    fit <- NA
    try(
      fit <- bbmle::mle2(flog, list(.a = as.numeric(skies[i,1]),
                                    .b = as.numeric(skies[i,2]),
                                    .c = as.numeric(skies[i,3]),
                                    .d = as.numeric(skies[i,4]),
                                    .e = as.numeric(skies[i,5])), method = method),
      silent = TRUE
    )

    if (any(try(fit@details$convergence, silent = TRUE), is.na(fit))) {
      fit <- NA
      coef <- as.numeric(skies[i,1:5])
      names(coef) <- c(".a", ".b", ".c", ".d", ".e")
      pred <- .cie_sky_model(AzP, Zp, AzS, Zs,
                             .a = as.numeric(skies[i,1]),
                             .b = as.numeric(skies[i,2]),
                             .c = as.numeric(skies[i,3]),
                             .d = as.numeric(skies[i,4]),
                             .e = as.numeric(skies[i,5]))
    } else {
      coef <- fit@coef
      coef[c(3,5)] <- abs(coef[c(3,5)])
      coef[c(2,4)] <- -abs(coef[c(2,4)])
      pred <- .cie_sky_model(AzP, Zp, AzS, Zs,
                             .a = fit@coef[1],
                             .b = fit@coef[2],
                             .c = fit@coef[3],
                             .d = fit@coef[4],
                             .e = fit@coef[5])
    }

    if (is.vector(custom_sky_coef)) {
      .start <- custom_sky_coef
    } else {
      .start <- as.numeric(skies[i,1:5]) %>% as.vector()
    }

    list(mle2_output = fit,
         coef = coef,
         obs = rl$sky_points$rl,
         pred = pred,
         zenith_dn = rl$zenith_dn,
         sun_coord = sun_coord,
         method = method,
         start = .start)
  }

  if (!is.null(std_sky_no)) {
    skies <- skies[std_sky_no,]
  }
  if (!is.null(general_sky_type) & is.null(std_sky_no)) {
    indices <- skies$general_sky_type == general_sky_type
    skies <- skies[indices,]
  }

  fit <- suppressWarnings(Map(.fun, 1:nrow(skies)))

  if (sun_coord$zenith_azimuth[1] > twilight) {
    indices <- match(7:15, as.numeric(rownames(skies)))
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

  .get_loglik <- function(x) {
    status <- if (length(x$coef) != 5) {
      "invalid_coef"
    } else if (!inherits(x$mle2_output, "mle2")) {
      "no_fit"
    } else if (x$mle2_output@details$convergence != 0) {
      "no_convergence"
    } else {
      "valid"
    }

    switch(status,
           invalid_coef = 10^10,
           no_fit = 10^9,
           no_convergence = 10^8,
           valid = x$mle2_output %>% summary() %>% .@m2logL %>%
             suppressWarnings())
  }

  error <- Map(.get_loglik, fit) %>% unlist()
  i <- which.min(error)
  if (error[i] > 10^7) {
    error <- Map(function(x) .calc_rmse(x$pred - x$obs), fit) %>% unlist()
  }

  model <- fit[[i]]
  if (model$sun_coord$zenith_azimuth[1] >= 90) {
      model$sun_coord$row_col <- c(NA, NA)
  }
  model
}

