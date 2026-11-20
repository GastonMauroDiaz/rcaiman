#' Distance to black pixels
#'
#' @description
#' Compute the Chebyshev distance from each `TRUE` pixel to the nearest
#' `FALSE` pixel in a binary image (often depicted in black).
#'
#' @details
#' Pixels marked `FALSE` are assigned distance `0`. For pixels marked `TRUE`,
#' the returned value is the Chebyshev distance to the nearest `FALSE` pixel,
#' using 8-neighbour connectivity. Larger values indicate pixels located deeper
#' inside connected `TRUE` regions and farther from region boundaries.
#'
#' If `bin` contains no `FALSE` pixels, distances are undefined. In this case,
#' a warning is issued and a raster filled with `NA` values is returned.
#'
#' @param bin logical [terra::SpatRaster-class] of one layer.
#'
#' @return Integer [terra::SpatRaster-class] of one layer named
#'   `"distance_to_black"`. Pixels corresponding to `FALSE` values in `bin` are
#'   assigned `0`; pixels corresponding to `TRUE` values receive their distance
#'   to the nearest `FALSE` pixel.
#'
#' @section Sky sampling:
#' Threshold this image to create candidate masks for [optimize_sky_mask()].
#' For integer `k`, `distance > k` removes TRUE pixels within `k` pixels of
#' black, in the Chebyshev metric. This is equivalent to repeated erosion
#' with a 3 by 3 square neighborhood, without treating the image exterior as
#' black. Only FALSE pixels within the image act as distance sources.
#'
#' After choosing a mask, its distances can also prioritize candidates in
#' [optimize_sampling_distance()] and [rem_nearby_indices()]: extract them
#' with [extract_dn()] and pass the resulting column with `prioritize = "higher"`.
#' Use the same priority for distance estimation and final selection.
#'
#' @seealso [label_components()], [optimize_sky_mask()],
#'   [optimize_sampling_distance()], [rem_nearby_indices()]
#'
#' @export
#'
#' @examples
#' bin <- set_rcaiman_geometry(terra::rast(nrows = 12, ncols = 16))
#' bin[] <- FALSE
#' bin[3:10, 3:14] <- TRUE
#' distance <- distance_to_black(bin)
#' plot(distance)
#'
#' # Progressively shrink the TRUE region, preserving the original at k = 0.
#' bin_list <- lapply(0:3, function(k) distance > k)
#' plot(terra::rast(bin_list))
#'
#' # Attach the distance as a priority value for candidate sampling indices.
#' sampling_indices <- expand.grid(col = 1:16, row = 1:12)
#' sampling_indices <- rem_unselected(sampling_indices, bin)
#' sampling_indices <- extract_dn(distance, sampling_indices, use_window = FALSE)
#' retained <- rem_nearby_indices(sampling_indices, min_dist = 3,
#'   col_id = "distance_to_black", prioritize = "higher", space = "planar")
#' # See optimize_sky_mask() for selection with angular separation.
distance_to_black <- function(bin) {
  .assert_logical_mask(bin)

  x <- terra::as.matrix(bin, wide = TRUE)

  distance <- distance_to_black_cpp(x)

  distance <- terra::setValues(bin, distance)
  names(distance) <- "distance_to_black"
  distance
}
