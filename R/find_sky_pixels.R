#' Find sky pixles
#'
#' Find sky pixels automatically
#'
#' This function assumes that (1) there is at least one pure sky pixel at the
#' level of cells of \eqn{30 \times 30} degrees, and (2) sky pixels have a
#' digital number (DN) greater than canopy pixels have.
#'
#' For each cell, it compute a quantile value and use it as a threshold to
#' select the pure sky pixels of the cell, which produce binarized image as a
#' result in a regional binarization fashion
#' (\code{\link{regional_thresholding}}). This process start with a quantile
#' probability of 0.99. After producing the binarized image, this function use a
#' search grid with cells of \eqn{5 \times 5} degrees to count how many cells on
#' the binarired image have at least one sky pixel. If the count does not reach
#' argument \code{no_of_samples}, it goes back to the binarization step but
#' decreasing the probability by 0.01 points.
#'
#' @inheritParams fit_coneshaped_model
#' @param no_of_samples Numeric vector of length one. Minimum number of samples
#'   required.
#'
#' @family MBLT functions
#'
#' @export
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}. This layer masks pixels that are very likely pure sky pixels.
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(caim$Blue)
#' bin <- find_sky_pixels(blue, z, a)
#' plot(bin)
#' }
find_sky_pixels <- function(r, z, a, no_of_samples = 30) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  stopifnot(length(no_of_samples) == 1)

  g30 <- sky_grid_segmentation(z, a, 30)
  g5 <- sky_grid_segmentation(z, a, 5)

  m <- mask_hs(z, 80, 90)

  prob <- 1
  count <- 0
  while (count <= no_of_samples) {
    if (prob < 0.9) {
      prob <- 1
      no_of_samples <- round(no_of_samples * 0.9)
      if (no_of_samples < 30) {
        stop(paste(
          "The function is not working properly.",
          "The problem might be related to inputs.",
          "Please, make sure they are OK."
        ))
      }
    }
    prob <- prob - 0.025
    bin <- regional_thresholding(r, g30, "Diaz2018", 0, 1, prob)
    bin[m] <- 0
    max_per_cell <- extract_feature(bin, g5, max, return_raster = FALSE)
    count <- sum(max_per_cell)
  }

  prob <- prob + 0.025
  while (count <= no_of_samples) {
    if (prob < 0.895) {
      stop("please, report error 001")
    }
    prob <- prob - 0.01
    bin <- regional_thresholding(r, g30, "Diaz2018", 0, 1, prob)
    bin[m] <- 0
    max_per_cell <- extract_feature(bin, g5, max, return_raster = FALSE)
    count <- sum(max_per_cell)
  }

  bin[is.na(z)] <- 0
  as.logical(bin)
}
