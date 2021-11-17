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
#' This pipeline is based on these two studies:
#' \insertCite{Lang2010;textual}{rcaiman} and
#' \insertCite{Diaz2018;textual}{rcaiman}.
#'
#' Details about how the original method by
#' \insertCite{Diaz2018;textual}{rcaiman} was adapted can be read here
#' \code{\link{ootb_mblt}}.
#'
#' The main differences between the original method by
#' \insertCite{Lang2010;textual}{rcaiman} and the one implemented here are: (1)
#' it is fully automatic, (2) the residuals of the CIE sky model
#' ($residuals=model-data$) are interpolated instead of the sky digital numbers
#' (the data), and (3) the final sky reconstruction is obtained by subtracting
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
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' sky <- ootb_sky_reconstruction(blue, z, a)
#' plot(sky)
#' ratio <- blue / sky
#' plot(ratio)
#' hist(ratio)
#' }
ootb_sky_reconstruction <- function(r, z, a) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  bin <- ootb_mblt(r, z, a)$bin
  bin[mask_hs(z, 80, 90)] <- 0
  g <- sky_grid_segmentation(z, a, 10)
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
  residu_i <- interpolate_dns(residu, sky_marks,
                              k = 3,
                              p = 2,
                              rmax = ncol(r) / 7,
                              use_window = TRUE)

  .sand_residu_i <- function(r, sky_cie, residu_i) {
    sky <- sky_cie - residu_i
    sky <- cover(sky, sky_cie)
    ratio <- r / sky
    current <- max(ratio[], na.rm = TRUE)
    rot <- 0 # rate of change
    no_loops <- 0
    while (current > 1.5 & rot <= 0 & no_loops < 30) {
      previous <- max(ratio[], na.rm = TRUE)
      no_loops <- no_loops + 1
      res_range <- range(residu_i[], na.rm = TRUE)
      if (rot == 0) index <- which.max(abs(res_range))
      residu_i[residu_i == res_range[index]] <- res_range[index] * 0.9
      sky <- sky_cie - residu_i
      sky <- cover(sky, sky_cie)
      ratio <- r / sky
      current <- max(ratio[], na.rm = TRUE)
      rot <- (current / previous) - 1
    }
    residu_i
  }
  delta <- abs(residu_i - .sand_residu_i(r, sky_cie, residu_i))
  delta_range <- range(delta[], na.rm = TRUE)
  if (sd(delta_range) != 0) {
    delta_mask <- delta > mean(delta_range)
    indices <- delta_mask[cellFromRowCol(r, sky_marks$row, sky_marks$col)]
    sky_marks <- sky_marks[!indices, ]
    model <- fit_cie_sky_model(r, z, a, sky_marks, sun_mark,
                               std_sky_no = NULL,
                               general_sky_type = NULL ,
                               use_window = TRUE,
                               twilight = TRUE,
                               rmse = FALSE,
                               method = "BFGS")
    sky_cie <- model$relative_luminance * model$zenith_dn
    residu <- sky_cie - r
    residu_i <- interpolate_dns(residu, sky_marks,
                                k = 3,
                                p = 2,
                                rmax = ncol(r) / 7,
                                use_window = TRUE)
  }
  sky <- sky_cie - residu_i
  cover(sky, sky_cie)
}
