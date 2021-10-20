#' CIE sky model
#'
#' Written by Gaston Mauro Diaz based on Pascal code by Mait Lang.
#'
#' Angles should be provided in radians.
#'
#' @param AzP Numeric vector.
#' @param Zp Numeric vector.
#' @param AzS Numeric vector of length one.
#' @param Zs Numeric vector of length one.
#' @param .a,.b,.c,.d,.e Numeric vector of length one. Sky model parameter.
#'
#' @references http://dx.doi.org/10.1016/j.energy.2016.02.054
#'
#' @return Numeric vector of length equal to AzP length.
.cie_sky_model <- function(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e) {
  .fun <- function(i) {
    # calculate angular distance between sky point and Sun
    delta_Az <- abs(AzP[i] - AzS)
    Chi <- acos(cos(Zs) * cos(Zp[i]) + sin(Zs) * sin(Zp[i]) * cos(delta_Az))

    # Gradation function
    Phi_Z <- 1 + .a * exp(.b / cos(Zp[i]))
    Phi_0 <- 1 + .a * exp(.b)
    gradation <- Phi_Z / Phi_0

    # Indicatrix function
    F_Chi <-  1 + .c * (exp(.d * Chi) - exp(.d * pi/2)) + .e * cos(Chi)^2
    F_Zs  <-  1 + .c * (exp(.d *  Zs) - exp(.d * pi/2)) + .e * cos(Zs)^2
    indicatrix <- F_Chi / F_Zs

    unname(gradation * indicatrix)
  }
  Map(.fun(i), seq_along(AzP)) %>% unlist()
}


CIE_sky_model_raster <- function(z, a, sun_coord, sky_coef) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(length(sun_coord) == 2)
  stopifnot(length(sky_coef) == 5)

  Zp <- degree2radian(z[])
  AzP <- degree2radian(a[])

  Zs <- degree2radian(sun_coord[1])
  AzS <- degree2radian(sun_coord[2])

  relative_luminance <- .cie_sky_model(AzP, Zp, AzS, Zs,
                                       sky_coef[1],
                                       sky_coef[2],
                                       sky_coef[3],
                                       sky_coef[4],
                                       sky_coef[5])


  z[] <- relative_luminance
  z
}


fit_cie_sky_model <- function(z, a, sky_marks, sun_mark, method = "SANN") {
  if (!requireNamespace("bbmle", quietly = TRUE)) {
    stop(paste("Package \"bbmle\" needed for this function to work.",
               "Please install it."
    ),
    call. = FALSE)
  }

  sun_coord <- degree2radian(c(a[sun_mark[2], sun_mark[1]],
                               z[sun_mark[2], sun_mark[1]]))
  AzS <- sun_coord[1]
  Zs <- sun_coord[2]


  cells <- cellFromRowCol(a, sky_marks$row, sky_marks$col)
  sky_marks$a <- a[cells]
  sky_marks$z <- z[cells]

  AzP <- degree2radian(sky_marks$a)
  Zp <- degree2radian(sky_marks$z)

  z_thr <- 5
  zenith_dn <- c()
  while (length(zenith_dn) < 20) {
    zenith_dn <- sky_marks$dn[sky_marks$z < z_thr]
    z_thr <- z_thr + 2
  }
  if ((z_thr - 2) > 20) warning(paste0("Zenith DN were estimated from pure sky",
                                 "pixels from zenith angle up to", z_thr, "."))

  zenith_dn <- mean(zenith_dn)
  sky_marks$dn <- sky_marks$dn / zenith_dn


  flog <- function(.a, .b, .c, .d, .e, S) {
    media <- .cie_sky_model(AzP, Zp, AzS, Zs, .a, .b, .c, .d, .e)
    - sum(dnorm(sky_marks$dn, mean = media, sd = exp(S),log = FALSE))
  }


  mm <-  bbmle::mle2(flog, list(.a = 4,
                                .b = 0,
                                .c = 2,
                                .d = -1,
                                .e = 0,
                                S = 0), method = method)


}

