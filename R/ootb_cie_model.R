#' Out-of-the-box CIE sky model
#'
#' @description
#' Fit a CIE general sky model to data obtained from a canopy image with an
#' histogram-based method.
#'
#' @note
#' This function is part of a paper under preparation.
#'
#' @param r numeric [terra::SpatRaster-class] of one layer. Typically, the blue
#'   band of a a canopy photograph. Digital numbers should be linearly related
#'   to radiance. See [read_caim_raw()] for details.
#' @param method character vector. Optimization methods for [fit_cie_model()]
#'   and [optim_sun_angles()].
#'
#' @inheritParams compute_canopy_openness
#' @inheritParams skygrid_segmentation
#' @inheritParams fit_cie_model
#' @inheritParams apply_by_direction
#'
#'
#' @return List with the following components:
#' \describe{
#'   \item{`rr`}{The input `rr` with an added `pred` column in
#'     `sky_points`, containing predicted values.}
#'   \item{`opt_result`}{List returned by [stats::optim()].}
#'   \item{`coef`}{Numeric vector of length five. CIE model coefficients.}
#'   \item{`sun_angles`}{Numeric vector of length two. Sun zenith and azimuth
#'     (degrees).}
#'   \item{`method`}{Character vector of length one. Optimization method used.}
#'   \item{`start`}{Numeric vector of length five. Starting parameters.}
#'   \item{`metric`}{Numeric value. Mean squared deviation as in
#'     \insertCite{Gauch2003;textual}{rcaiman}.}
#'   \item{`method_sun`}{Character vector of length one or `NULL`.
#'     Method used to optimice sun angles.}
#' }
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' model <- ootb_cie_model(r, z, a, m)
#'
#' plot(model$rr$sky_points$rr, model$rr$sky_points$pred)
#' abline(0,1)
#' lm(model$rr$sky_points$pred~model$rr$sky_points$rr) %>% summary()
#'
#' sky <- cie_image(z, a, model$sun_angles, model$coef) * model$rr$zenith_dn
#' plot(sky)
#' ratio <- r/sky
#' plot(ratio)
#' plot(ratio > 1.05)
#' plot(ratio > 1.15)
#' }
ootb_cie_model <- function(r, z, a, m,
                           method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
                           parallel = TRUE,
                           cores = NULL,
                           leave_free = 1,
                           logical = TRUE
) {

  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }
  #more basic checks are handled by the functions called below

  sky <- apply_by_direction(r, z, a, m,
                            spacing = 10, laxity = 2.5,
                            fov = 40,
                            method = "detect_bg_dn",
                            parallel = parallel,
                            cores = cores)

  g <- skygrid_segmentation(z, a, 22.5, first_ring_different = FALSE)

  sky_points <- sample_sky_points(sky$n, !is.na(sky$dn), g,
                                   dist_to_black = 1)

  rr <- extract_rr(sky$dn, z, a, sky_points,
                   no_of_points = 3,
                   use_window = FALSE)

  # Sun coordinate
  sun_angles <- estimate_sun_angles(r, z, a, m, NULL, NULL,
                                    method = "assume_veiled")

  # Force the sun low and add that to the results
  civic_twilight <- c(seq(sun_angles["z"], 90, length = 5), seq(91, 96, 1))

  suns <- lapply(seq_along(civic_twilight),
                 function(i) c(z = civic_twilight[i], sun_angles["a"]))

  # Model
  if (parallel) {
    # Only to avoid note from check, code is OK without this line.
    sun <- NA

    .with_cluster(cores, {
      models <- foreach(sun = suns) %dopar% {
        fit_cie_model(rr, sun, method = method)
      }
    })
  } else {
    models <- lapply(suns, function(sun) fit_cie_model(rr, sun,
                                                       method = method))
  }

  metrics <- lapply(seq_along(models), function(i) models[[i]]$metric) %>%
    unlist
  i <- which.min(metrics)
  model <- models[[i]]

  model <- optim_sun_angles(model, method)
  model
}
