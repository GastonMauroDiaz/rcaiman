#' Remove small gaps
#'
#' @description
#' Remove small canopy gaps by applying a binary morphological opening with a
#' diamond-shaped structuring element. The operation can also remove narrow
#' parts and protrusions from larger gaps.
#'
#' @details
#' Opening is an erosion followed by a dilation with the same structuring
#' element. Both steps use [terra::focal()] on the pixel grid. The diamond
#' includes offsets whose absolute row and column displacements sum to at most
#' `(size - 1) / 2` pixels.
#'
#' At image edges, only pixels within the image participate. Outside positions
#' are neutral: TRUE for erosion and FALSE for dilation. The output retains
#' the original grid, including when the structuring element exceeds the image
#' dimensions. With `size = 1`, the input is returned unchanged.
#'
#' For selection of whole components by size or shape, use [label_components()]
#' and the descriptors provided by [compute_features()].
#'
#' @param bin logical, NA-free [terra::SpatRaster-class] of one layer.
#' @param size positive odd integer of length one. Width and height of the
#'   diamond's bounding square, in pixels.
#'
#' @return Logical, NA-free [terra::SpatRaster-class] on the same grid as `bin`.
#' @seealso [distance_to_black()], [rem_isolated_black_pixels()]
#' @export
#' @examples
#' bin <- set_rcaiman_geometry(terra::rast(nrows = 12, ncols = 16))
#' bin[] <- FALSE
#' bin[3:10, 7:14] <- TRUE
#' bin[3, 3] <- TRUE
#' cleaned <- rem_small_gaps(bin, size = 3)
#' plot(c(bin, cleaned))
rem_small_gaps <- function(bin, size = 9) {
  .assert_logical_mask(bin)
  .check_vector(size, "odd_integerish", 1, sign = "positive")
  if (size == 1) return(bin)

  radius <- (size - 1) / 2
  offsets <- abs(seq.int(-radius, radius))
  kernel <- outer(offsets, offsets, "+")
  kernel <- ifelse(kernel <= radius, 1, NA_real_)

  step <- function(x, fun, outside) {
    # Padding also supports rasters smaller than the focal window.
    padded <- terra::extend(x, c(radius, radius), fill = outside)
    filtered <- terra::focal(padded, w = kernel, fun = fun,
                             na.rm = TRUE, fillvalue = outside)
    terra::crop(filtered, bin)
  }
  eroded <- step(bin, "min", 1)
  opened <- step(eroded, "max", 0)
  opened <- as.logical(opened)
  names(opened) <- names(bin)
  opened
}
