#' Find sky pixels
#'
#' Find sky pixels automatically
#'
#' This function assumes that (1) there is at least one pure sky pixel at the
#' level of cells of \eqn{30 \times 30} degrees, and (2) sky pixels have a
#' digital number (DN) greater than canopy pixels have.
#'
#' For each \eqn{30 \times 30} cell, this method computes a quantile value and
#' use it as a threshold to select the pure sky pixels of the given cell. As a
#' result, a binarized image is produced in a regional binarization fashion
#' (\code{\link{regional_thresholding}}). This process start with a quantile
#' probability of 0.99. After producing the binarized image, this function use a
#' search grid with cells of \eqn{5 \times 5} degrees to count in how many of
#' these cells are at least one sky pixel (pixels equal to one in the binarized
#' image). If the percentage of  cells with sky pixels does not reach argument
#' \code{sample_size_pct}, it goes back to the binarization step but decreasing
#' the probability by 0.01 points.
#'
#' If probability reach 0.9 and the \code{sample_size_pct} criterion were not
#' yet satisfied, the \code{sample_size_pct} is decreased one percent and the
#' process start all over again.
#'
#'
#'
#'
#' @inheritParams ootb_mblt
#' @param sample_size_pct Numeric vector of length one. Minimum sample size
#'   percentage required. The population is comprised of 1296 cells of \eqn{5
#'   \times 5} degrees.
#'
#' @family MBLT functions
#'
#' @export
#' @return Object from class list containing an object of class
#'   \linkS4class{SpatRaster} (named ‘bin’), and two numeric vector of length
#'   one (named ‘prob’ and ‘sample_size_pct’). Object ‘bin’ masks pixels that
#'   are very likely pure sky pixels. The objects ‘prob’ and ‘sample_size_pct’
#'   are described in the Details section.
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(caim$Blue)
#' bin <- find_sky_pixels(blue, z, a)$bin
#' plot(bin)
#' }
find_sky_pixels <- function(r, z, a, sample_size_pct = 30) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  stopifnot(length(sample_size_pct) == 1)

  if (sample_size_pct > 30) sample_size_pct <- 31

  cells_number <- (360/5) * (90/5)

  g30 <- sky_grid_segmentation(z, a, 30)
  g5 <- sky_grid_segmentation(z, a, 5)

  m <- mask_hs(z, 80, 90)

  prob <- 1
  count <- 0
  while ((count / cells_number) * 100 <= sample_size_pct) {
    if (prob < 0.9) {
      prob <- 1
      sample_size_pct <- sample_size_pct - 1
      if (sample_size_pct < 5) {
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
  while ((count / cells_number) * 100 <= sample_size_pct) {
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
  list(bin = bin, prob = prob, sample_size_pct = sample_size_pct)
}
