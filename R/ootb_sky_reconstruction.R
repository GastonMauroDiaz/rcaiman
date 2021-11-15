#' Out-of-the-box sky reconstruction
#'
#' Build an above canopy image from a single below canopy image.
#'
#' This function is a hard-coded version of a pipeline that uses these main
#' functions: \code{\link{ootb_mblt}}, \code{\link{fit_cie_sky_model}}, and
#' \code{\link{interpolate_dns}}. The code can be easily inspected by calling
#' \code{ootb_sky_reconstruction} --no parenthesis. Advanced users could use
#' that code as a template.
#'
#' This pipeline integrate these two studies:
#' \insertCite{Lang2010;textual}{rcaiman} and
#' \insertCite{Diaz2018;textual}{rcaiman}.
#'
#' Details about how the original method by
#' \insertCite{Diaz2018;textual}{rcaiman} was adapted can be read here
#' \code{\link{ootb_mblt}}.
#'
#' The main difference between the original method by
#' \insertCite{Lang2010;textual}{rcaiman} and the one implemented here are: (1)
#' it is fully automatic, (2) the interpolation is done in the raster space of a
#' cylindrical projection instead of the raster space of an hemispherical
#' equiangular projection, (3) the residuals of the CIE sky model
#' ($residuals=model-data$) are interpolated instead of the sky digital numbers
#' (the data), and (4) the final sky reconstruction is obtained by subtracting
#' the interpolated residuals to the CIE sky model instead of by calculating a
#' weighted average parameterized by the user.
#'
#' The recommended input for this function is data pre-processed with the HSP
#' software package \insertCite{Lang2013}{rcaiman}. Please, refer to
#' \code{\link{write_sky_marks}} for additional details about HSP, and refer to
#' \code{\link{fit_cie_sky_model}} and \code{\link{interpolate_dns}} to know
#' why the HSP pre-processing is convenient.
#'
#'
#' @inheritParams ootb_mblt
#' @inheritParams zenith_image
#'
#'
#' @export
#'
#' @references \insertRef{Lang2010}{rcaiman}
#'
#'   \insertRef{Lang2013}{rcaiman}
#'
#'   \insertRef{Diaz2018}{rcaiman}
#'
#'
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
#' blue <- gbc(r$Blue)
#' sky <- ootb_sky_reconstruction(blue, lens("Nikon_FCE9"))
#' plot(sky)
#' ratio <- blue / sky
#' plot(ratio)
#' hist(ratio)
#' plot(ratio > 1.3)
#' ratio[ratio > 1.3] <- 1.3
#' plot(ratio)
#' }
ootb_sky_reconstruction <- function(r, lens_coef) {
  stopifnot(ncol(r) == nrow(r))
  .check_if_r_was_normalized(r)
  z <- zenith_image(ncol(r), lens_coef)
  if (.get_max(z) > 90) stop(paste("Please check your \"lens_coef\" input",
                                   "with \"test_lens_coef()\"."))
  a <- azimuth_image(z)
  g <- sky_grid_segmentation(z, a, 10)
  bin <- ootb_mblt(r, z, a)$bin

  m <- mask_hs(z, 80, 90)
  bin[m] <- 0

  sky_marks <- extract_sky_marks(r, bin, g,
                                 dist_to_plant = 3,
                                 min_raster_dist = 3)
  sun_mark <- extract_sun_mark(r, bin, z, a, g, max_angular_dist = 45)
  model <- fit_cie_sky_model(r, z, a, sky_marks, sun_mark,
                             std_sky_no = NULL,
                             general_sky_type = NULL ,
                             use_window = TRUE,
                             twilight = TRUE,
                             rmse = FALSE,
                             method = "BFGS")
  sky_cie <- model$relative_luminance * model$zenith_dn
  residu <- sky_cie - r
  residu_i <- interpolate_reproj(residu, z, a, model$coef, sky_marks,
                                 k = 3,
                                 p = 2,
                                 rmax = 20,
                                 use_window = TRUE)
  sky <- sky_cie - residu_i
  cover(sky, sky_cie)
}
