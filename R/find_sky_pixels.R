#' Find sky pixels
#'
#' Find sky pixels automatically.
#'
#' This function assumes that:
#'
#' * there is at least one pure sky pixel at the level of cells of \eqn{30
#' \times 30} degrees, and
#' * sky pixels have a digital number (DN) greater than canopy pixels have.
#'
#' For each \eqn{30 \times 30} cell, this method computes a quantile value and
#' uses it as a threshold to select the pure sky pixels from the given cell. As
#' a result, a binarized image is produced in a regional binarization fashion
#' ([regional_thresholding()]). This process starts with a quantile
#' probability of 0.99. After producing the binarized image, this function uses
#' a search grid with cells of \eqn{5 \times 5} degrees to count how many of
#' these cells have at least one sky pixel (pixels equal to one in the binarized
#' image). If the percentage of  cells with sky pixels does not reach argument
#' `sample_size_pct`, it goes back to the binarization step but decreasing
#' the probability by 0.01 points.
#'
#' If probability reach 0.9 and the `sample_size_pct` criterion were not
#' yet satisfied, the `sample_size_pct` is decreased one percent and the
#' process starts all over again.
#'
#' @inheritParams ootb_mblt
#' @param sample_size_pct Numeric vector of length one. Minimum percentage of
#'   cells to sample. The population is comprised of 1296 cells of \eqn{5 \times
#'   5} degrees.
#'
#' @family Tool Functions
#'
#' @export
#' @return An object of class [SpatRaster-class] with values `0` and
#'   `1`. This layer masks pixels that are very likely pure sky pixels.
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize(., 0, 20847)
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' r <- caim$Blue
#' r[is.na(r)] <- 0
#' bin <- find_sky_pixels(r, z, a)
#' plot(bin)
#' }
find_sky_pixels <- function(r, z, a, sample_size_pct = 30) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  .was_normalized(r)
  stopifnot(length(sample_size_pct) == 1)

  if (sample_size_pct > 60) {
    warning("A large \"sample_size_pct\" may demand long processing time")
  }

  cells_number <- (360/5) * (90/5)

  g30 <- sky_grid_segmentation(z, a, 30)
  g5 <- sky_grid_segmentation(z, a, 5)

  m <- select_sky_vault_region(z, 80, 90)

  prob <- 1
  count <- 0
  while ((count / cells_number) * 100 <= sample_size_pct) {
    if (prob < 0.9) {
      prob <- 1
      warning(paste0("'sample_size_pct' was forced to ", sample_size_pct))
      sample_size_pct <- sample_size_pct - 1
      if (sample_size_pct < 5) {
        stop("The 'sample_size_pct' can not be forced below 5")
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
      stop("'prob' reached 0.895")
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
