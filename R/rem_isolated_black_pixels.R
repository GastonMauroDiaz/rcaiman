#' Remove isolated black pixels
#'
#' @description
#' Replace single black pixels (`FALSE`) that are fully surrounded by white
#' pixels (`TRUE`) with white. Uses 8-connectivity.
#'
#' @inheritParams write_bin
#'
#' @returns Logical [terra::SpatRaster-class] of one layer.
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' r <- caim$Blue
#'
#' bin <- binarize_by_region(r, ring_segmentation(z, 15), "thr_isodata")
#'
#' plot(bin)
#' plot(rem_isolated_black_pixels(bin))
#' }
rem_isolated_black_pixels <- function(bin) {
  .assert_logical_mask(bin)

  # Define a 3x3 Laplacian-like kernel
  ma <- matrix(c(1, 1, 1,
                 1, -8, 1,
                 1, 1, 1), ncol = 3)

  # Identify isolated black pixels: they have 8 white neighbors
  isolated <- terra::focal(bin, w = ma) == 8

  # Set them to white
  bin[isolated] <- 1

  # Ensure return values are logical (TRUE/FALSE)
  as.logical(bin)
}
