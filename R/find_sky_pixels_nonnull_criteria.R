#' Find sky pixels following the non-null criteria
#'
#' This function assumes that (1) there is at least one pure sky pixel at the
#' level of cells of \eqn{30 \times 30} degrees, and (2) sky pixels have a
#' digital number (DN) greater than canopy pixels have.
#'
#' A sky grid of \eqn{30 \times 30} degrees along with the arguments \code{r},
#' \code{prob} and \code{slope} are passed to
#' \code{\link{regional_thresholding}} in order to produce a binarized image.
#' Then, a sky grid of \eqn{10 \times 10} degrees is used to compute the number
#' of cells having none sky pixels (the so-called null cells).
#'
#' The process is repeated but increasing \code{slope} in steps of 0.05 till the
#' the number of null cells increases.
#'
#' The \code{prob} argument can be obtained with \code{\link{find_sky_pixels}}.
#'
#' @inheritParams find_sky_pixels
#' @inheritParams regional_thresholding
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(caim$Blue)
#' bin <- find_sky_pixels_nonnull_criteria(blue, z, a)
#' plot(bin)
#' }
find_sky_pixels_nonnull_criteria <- function(r, z, a, prob = 1, slope = 0.5) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  g30 <- sky_grid_segmentation(z, a, 30)
  g <- sky_grid_segmentation(z, a, 10)
  g[mask_hs(z, 0, 10) | mask_hs(z, 70, 90)] <- NA
  .get_no_of_nulls <- function(bin) {
    no_of_nulls <- extract_feature(bin, g, mean, return_raster = FALSE)
    sum(no_of_nulls == 0)
  }
  slp <- 0.5
  unlock <- TRUE
  while (unlock) {
    no_nulls_1 <- .get_no_of_nulls(regional_thresholding(r, g30, "Diaz2018",
                                                         0, slp, prob))
    slp <- slp + 0.05
    suppressWarnings(rm("no_nulls_2"))
    try(
      no_nulls_2 <- .get_no_of_nulls(regional_thresholding(r, g30, "Diaz2018",
                                                           0, slp, prob)),
      silent = TRUE
    )
    if (exists("no_nulls_2")) {
      unlock <- no_nulls_2 <= no_nulls_1
    } else {
      unlock <- FALSE
    }
    if (!unlock) slp <- slp - 0.05
  }
  regional_thresholding(r, g30, "Diaz2018", 0, slp, prob)

}
