#' Fit CIE sky model
#'
#' Use general-purpuse optimization to find the parameters of the CIE sky model
#' that best fit data sampled from a canopy photograph.
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
#' r <- read_caim("manipulate/IMG_1013.pgm") %>% normalize_minmax()
#' z <- zenith_image(ncol(r), lens())
#' a <- azimuth_image(z)
#' manual_input <- read_manual_input(".", "IMG_1013" )
#' sun_coord <- manual_input$sun_coord$row_col
#' sun_coord <- zenith_azimuth_from_row_col(z, sun_coord, lens())
#' sky_points <- manual_input$sky_points
#' rr <- extract_rel_radiance(r, z, a, sky_points)
#' model <- fit_cie_sky_model(rr, sun_coord)
#' cie_sky <- model$relative_luminance * model$zenith_dn
#' plot(r/cie_sky)
#'
#' r <- read_caim("manipulate/IMG_1013.pgm")
#' sky_coef <- read_opt_sky_coef(".", "IMG_1013")
#' cie_sky_m <- cie_sky_image(z, a, sun_coord$zenith_azimuth, sky_coef)
#' cie_sky_m <- cie_sky_manual * manual_input$zenith_dn
#' plot(r/cie_sky_m)
#' ````
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
#' @param rr An object of class *list*. The output of [extract_rel_radiance()] or an
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
#'   these: `"Overcast"`, `"Clear"`, or `"Partly cloudy"`. See Table 1 from
#'   \insertCite{Li2016;textual}{rcaiman} for additional details.
#' @param twilight Numeric vector of length one. Sun zenith angle (in degrees).
#'   If the sun zenith angle provided through the `sun_coord` argument is below
#'   this value, sun zenith angles from 85 to 96 degrees will be tested when
#'   selecting initial optimization parameters from the general sky types
#'   _Clear_ or _Partially Cloudy_ (specifically, for standard sky numbers 7 to
#'   15). This adjustment is necessary because [extract_sun_coord()] can
#'   mistakenly identify the visible center of the solar corona as the solar
#'   disk. Since [extract_sun_coord()] cannot output a zenith angle below 90
#'   degrees, setting this value to 90 is equivalent to disabling this step.
#'   Actually, [extract_sun_coord()] cannot output a value very close to 90,
#'   therefore the testing start at 85, civic twilight is from 90 to 96 degrees.
#' @inheritParams stats::optim
#' @param loss Character vector of length one. Specifies the error metric to
#'   use. Options are `"MAE"` (Mean Absolute Error) or `"RMSE"` (Root Mean
#'   Squared Error). Defaults to `"MAE"`. `"MAE"` is more robust to inaccurate
#'   sky points and outliers, such as those found in cloudy skies.
#'
#' @references \insertAllCited{}
#'
#' @return A _list_ with the following components:
#' \itemize{
#'   \item The output produced by [stats::optim]
#'   \item The 5 coefficients of the CIE model
#'   \item Observed values
#'   \item Predicted values
#'   \item The estimated digital number at the zenith
#'   \item The sun coordinates
#'   \item The method used for optimization
#'   \item The starting parameters
#' }
#'
#' @family  Sky Reconstruction Functions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
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
#' rr <- extract_rel_radiance(caim$Blue, z, a, sky_points)
#'
#' set.seed(7)
#' model <- fit_cie_sky_model(rr, sun_coord,
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
fit_cie_sky_model <- function(rr, sun_coord,
                              custom_sky_coef = NULL,
                              std_sky_no = NULL,
                              general_sky_type = NULL ,
                              twilight = 60,
                              method = "BFGS",
                              loss = "MAE") {

  stopifnot(general_sky_type == "Overcast" |
            general_sky_type == "Partly cloudy" |
            general_sky_type == "Clear" |
            is.null(general_sky_type))
  stopifnot(is.data.frame(rr$sky_points))
  stopifnot(length(sun_coord$zenith_azimuth) == 2)

  # Manage the set of start parameter according to user choice
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
  if (!is.null(std_sky_no)) {
    skies <- skies[std_sky_no,]
  } else if (!is.null(general_sky_type)) {
    indices <- skies$general_sky_type == general_sky_type
    skies <- skies[indices,]
  }

  # Retrieve coordinates and convert to radians
  AzP <- .degree2radian(rr$sky_points$a)
  Zp <- .degree2radian(rr$sky_points$z)

  # Try all start parameters (brute force approach)
  .fun <- function(i) {
    # This has to be inside the function because of the twilight code chunk
    sun_a_z <- .degree2radian(rev(sun_coord$zenith_azimuth))
    AzS <- sun_a_z[1]
    Zs <- sun_a_z[2]

    # **Note 1**: when the .c or .e parameters are negative, the sun can be
    # darker than the sky, which does not make sense. So, the code avoids that
    # possibility by using abs()
    # **Note 2**: Positive values of b create unrealistic values at the
    # horizon. Positive values of d produce negative values, which are
    # unrealistc
    fcost <- function(params) {
      .a <- params[1]
      .b <- -abs(params[2])
      .c <- abs(params[3])
      .d <- -abs(params[4])
      .e <- abs(params[5])

      x <- .cie_sky_model(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e)
      residuals <- x - rr$sky_points$rr

      if (loss == "MAE") {
        return(stats::median(abs(residuals)))
      } else if (loss == "RMSE") {
        return(sqrt(mean(residuals^2)))
      } else {
        stop("The argument 'loss' must be 'MAE' or 'RMSE'")
      }
    }

    start_params <- skies[i, 1:5] %>% as.numeric()
    fit <- tryCatch(
      stats::optim(par = start_params,
            fn = fcost,
            method = method),
      error = function(e) list(par = start_params, convergence = NULL)
    )

    coef <- fit$par
    names(coef) <- c("a", "b", "c", "d", "e")
    coef[c(3,5)] <- abs(coef[c(3,5)])
    coef[c(2,4)] <- -abs(coef[c(2,4)])
    pred <- .cie_sky_model(AzP, Zp, AzS, Zs,
                           .a = coef[1],
                           .b = coef[2],
                           .c = coef[3],
                           .d = coef[4],
                           .e = coef[5])

    list(opt_result = fit,
         coef = coef,
         obs = rr$sky_points$rr,
         pred = pred,
         zenith_dn = rr$zenith_dn,
         sun_coord = sun_coord,
         method = method,
         start = start_params)
  }

  opt_result <- suppressWarnings(Map(.fun, 1:nrow(skies)))

  # Force the sun low and add that to the results
  civic_twilight <-  c(seq(85, 96, 1)) #a bit above horizon because bbox method
  if (sun_coord$zenith_azimuth[1] > twilight) {
    indices <- match(7:15, as.numeric(rownames(skies)))
    indices <- indices[!is.na(indices)]
    if (length(indices) != 0) {
      skies <- skies[indices,]
      for (i in seq_along(civic_twilight)) {
        sun_coord$zenith_azimuth <-  c(civic_twilight[i],
                                      sun_coord$zenith_azimuth[2])
        opt_result <- c(opt_result, suppressWarnings(Map(.fun, 1:nrow(skies))))
      }
    }
  }

  # Choose the best result
  .get_metric <- function(x) {
    residuals <- x$pred - x$obs
    if (loss == "MAE") {
      return(stats::median(abs(residuals)))
    } else if (loss == "RMSE") {
      return(sqrt(mean(residuals^2)))
    }
  }

  metric <- Map(.get_metric, opt_result) %>% unlist()
  i <- which.min(metric)
  model <- opt_result[[i]]

  if (model$sun_coord$zenith_azimuth[1] >= min(civic_twilight)) {
      model$sun_coord$row_col <- c(NA, NA)
  }
  model
}

