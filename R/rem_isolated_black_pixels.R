#' Remove isolated black pixels
#'
#' @description
#' Replace single black pixels (`FALSE`) that are fully surrounded by white
#' pixels (`TRUE`). Uses 8-connectivity.
#'
#' @inheritParams write_bin
#'
#' @returns Logical [terra::SpatRaster-class] of one layer.
#' @export
#'
#' @examples
#' bin <- matrix(TRUE, nrow = 8, ncol = 10)
#' bin[3, 3] <- FALSE           # One isolated black pixel.
#' bin[3, 7] <- FALSE           # Two black pixels touching diagonally.
#' bin[4, 8] <- FALSE
#' bin[6:7, 3:4] <- FALSE       # A larger black patch.
#' bin <- set_rcaiman_geometry(rast(bin)) > 0
#'
#' cleaned <- rem_isolated_black_pixels(bin)
#' plot(c(bin, cleaned), col = c("black", "white"),
#'      main = c("Original", "Isolated pixel removed"))
#'
#' # Only the isolated pixel changes; the diagonal pair and larger patch remain.
#' sum(cleaned[] != bin[])
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
