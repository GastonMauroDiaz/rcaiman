#' Defuzzify fuzzy classification
#'
#' This function translates degree of membership into Boolean logic using a
#' regional approach. The result will ensure that the fuzzy and Boolean version
#' will agree at the chosen level of aggregation (controlled by the argument
#' `segmentation`). This method makes perfect sense to translate a subpixel
#' classification of gap fraction (or a linear ratio) into a binary product.
#'
#' @note
#'
#' This method is also available in the HSP software package
#' \insertCite{Lang2013}{rcaiman}.
#'
#' @param mem An object of the class [SpatRaster-class]. Degree of
#'   membership.
#' @param segmentation An object of the class [SpatRaster-class], such as
#'   the result of a call to [sky_grid_segmentation()].
#'
#' @return An object of the class [SpatRaster-class] containing binary
#'   information.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' r <- correct_vignetting(r, z, c(0.0638, -0.101)) %>% normalize()
#' bin <- find_sky_pixels(r, z, a)
#' bin <- ootb_mblt(r, z, a, bin)
#' plot(bin$bin)
#' ratio <- r / bin$sky_s
#' ratio <- normalize(ratio, 0, 1, TRUE)
#' plot(ratio)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin2 <- defuzzify(ratio, g)
#' plot(bin2)
#' plot(abs(bin$bin - bin2))
#' }
defuzzify <- function (mem, segmentation) {
  .is_single_layer_raster(mem)
  .was_normalized(mem)
  mem[is.na(mem)] <- 0
  .is_single_layer_raster(segmentation)

  .fun <- function(x) {
    no_of_pixels <- round(mean(x) * length(x))
    if (no_of_pixels > 0 &
        no_of_pixels != length(x)) {
      indices <- order(x, decreasing = TRUE)[1:no_of_pixels]
      x[indices] <- 1
      x[x != 1] <- 0
    } else {
      x <- as.numeric(x > 0.5)
    }
    x
  }

  cells <- mem
  terra::values(cells) <- 1:ncell(mem)
  cells <- tapply(terra::values(cells),
                   terra::values(segmentation), function(x) x)
  cells <- unlist(cells)
  bin <- tapply(terra::values(mem), terra::values(segmentation), .fun)
  mem[cells] <- unlist(bin)
  as.logical(mem)
}
