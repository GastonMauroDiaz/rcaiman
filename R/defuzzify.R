#' Defuzzify a fuzzy classification
#'
#' This function translates degree of membership into Boolean logic using a
#' regional approach. The output ensures that the fuzzy and Boolean versions
#' remain consistent at the specified level of aggregation (controlled by the
#' argument `segmentation`). This method makes perfect sense to translate a
#' subpixel classification of gap fraction (or a linear ratio) into a binary
#' product.
#'
#' @note
#'
#' This method is also available in the HSP software package
#' \insertCite{Lang2013}{rcaiman}.
#'
#' @param mem A [SpatRaster-class] object representing the degree of membership
#'   in a fuzzy classification.
#' @param segmentation An object of the class [SpatRaster-class] such as an
#'   object returned by [sky_grid_segmentation()].
#'
#' @return An object of the class [SpatRaster-class] containing binary
#'   information.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path), z, a)
#'
#' sky <- ootb_build_sky_vault(r, z, a, ootb_sky$sky_points, ootb_sky)
#'
#' ratio <- r / sky$sky
#' ratio <- normalize_minmax(ratio, 0, 1, TRUE)
#' plot(ratio)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin2 <- defuzzify(ratio, g)
#' plot(bin2) # unsatisfactory results due to light conditions
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
