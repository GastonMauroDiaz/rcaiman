#' Optimize sun angular coordinates
#'
#' @inheritParams fit_cie_sky_model
#' @inheritParams cie_sky_image
#'
#' @return Numeric vector of length two, where the first element is the
#'   solar zenith angle and the second is the solar azimuth angle, both
#'   expressed in degrees.
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # See fit_cie_sky_model() for details on below file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' # x11()
#' # plot(caim$Blue)
#' # sun_zenith_azimuth <- click(c(z, a), 1) %>% as.numeric()
#' sun_zenith_azimuth <- c(49.5, 27.42481) #taken with above lines then hardcoded
#'
#' sun_row_col <- row_col_from_zenith_azimuth(z, a,
#'                                            sun_zenith_azimuth[1],
#'                                            sun_zenith_azimuth[2])
#' points(sun_row_col[2], nrow(caim) - sun_row_col[1], col = 3, pch = 1)
#'
#' rr <- extract_rel_radiance(caim$Blue, z, a, sky_points)
#'
#' set.seed(7)
#' model <- fit_cie_sky_model(rr, sun_zenith_azimuth,
#'                            general_sky_type = "Clear",
#'                            twilight = 90,
#'                            method = "BFGS")
#'
#' sun_zenith_azimuth <- optim_sun_zenith_azimuth(c(49.5, 27.42481), rr,
#'                                                model$coef, method = "CG")
#' sun_zenith_azimuth
#' sun_row_col <- row_col_from_zenith_azimuth(z, a,
#'                                            sun_zenith_azimuth[1],
#'                                            sun_zenith_azimuth[2])
#' points(sun_row_col[2], nrow(caim) - sun_row_col[1], col = 4, pch = 1)
#' }
optim_sun_zenith_azimuth <- function(sun_zenith_azimuth,
                                     rr,
                                     sky_coef,
                                     method = "BFGS") {
  .normalize_angles <- function(zenith, azimuth) {
    if (zenith < 0) {
      zenith <- -zenith        # Make zenith positive
      azimuth <- azimuth + 180 # Rotate azimuth 180°
    }

    # Ensure azimuth is in [0, 360)
    azimuth <- azimuth %% 360

    return(c(zenith, azimuth))
  }
  .refine_sun_coord <- function(param) {
    zenith <- param[1] * 9
    azimuth <- param[2] * 36

    chi <- calc_spherical_distance(zenith %>% .degree2radian(),
                                   azimuth %>% .degree2radian(),
                                   sun_zenith_azimuth[1] %>% .degree2radian(),
                                   sun_zenith_azimuth[2] %>% .degree2radian())


    pred <- .cie_sky_model(AzP = rr$sky_points$a %>% .degree2radian(),
                           Zp = rr$sky_points$z %>% .degree2radian(),
                           AzS = azimuth %>% .degree2radian(),
                           Zs =  zenith %>% .degree2radian(),
                           sky_coef[1], sky_coef[2], sky_coef[3],
                           sky_coef[4], sky_coef[5])

    (.calc_rmse(pred - rr$sky_points$rr) / mean(rr$sky_points$rr)) *
      max(1, .radian2degree(chi)/10)
  }
  fit <- tryCatch(stats::optim(c(sun_zenith_azimuth[1]/9,
                                 sun_zenith_azimuth[2]/36),
                               .refine_sun_coord,
                               method = method),
                  error = function(e) c(NA, NA))
  tryCatch(.normalize_angles(fit$par[1]*9, fit$par[2]*36),
    error = function(e) c(NA, NA))
}

