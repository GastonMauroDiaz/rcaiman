#' Defuzzify a fuzzy classification
#'
#' @description
#' Convert fuzzy membership values to a binary classification within segments,
#' approximating each segment's total membership by an integer count of pixels.
#'
#' @param mem numeric [terra::SpatRaster-class] of one layer. Membership values
#'   in `[0, 1]`. Missing values are treated as zero and remain `FALSE` in the
#'   output. Use a mask to exclude missing measurements from subsequent analyses.
#' @inheritParams compute_mean
#'
#' @details
#' Within each included segment, the function selects the pixels with the largest
#' memberships until their count equals `round(sum(mem))`. R's rounding rule
#' applies, including rounding exact halves to the nearest even integer. Ties
#' are resolved in raster cell order, from left to right and top to bottom.
#' The binary fraction therefore approximates the mean membership; an exact
#' match is not always possible with an integer number of pixels.
#'
#' This is an unweighted pixel count, not an area-weighted fraction. Missing
#' memberships contribute zero to the sum and are never selected. Values outside
#' `[0, 1]` produce a warning; the requested count is bounded by zero and the
#' number of nonmissing pixels in the segment. Membership values are not clipped.
#'
#' @return Logical [terra::SpatRaster-class] without `NA`, on the same grid as
#'   `mem`. Excluded labels (missing labels and, by default, zero) retain the
#'   initial value `FALSE`. This is a fill value, not a classification outside
#'   the analysis domain; delimit that domain with a mask.
#'
#' @note This method is also available in the HSP software package. See
#'   [hsp_compat()].
#' @seealso [binarize_by_region()]
#' @export
#' @examples
#' mem <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 4))
#' terra::values(mem) <- c(0.9, 0.6, 0.2, 1, 0.8, 0.3, NA, 1)
#' seg <- terra::rast(mem)
#' terra::values(seg) <- c(1, 1, 1, 0, 2, 2, 2, NA)
#' bin <- defuzzify(mem, seg)
#' plot(bin)
#' # Include label zero as another region when appropriate.
#' defuzzify(mem, seg, ignore_label_0 = FALSE)
defuzzify <- function(mem, seg, ignore_label_0 = TRUE) {
  .assert_single_layer(mem)
  .assert_same_geom(mem, seg)
  labels <- .feature_labels(seg, ignore_label_0)
  membership <- terra::values(mem, mat = FALSE)
  included <- !is.na(labels$index)
  .check_vector(membership[included], "numeric", allow_na = TRUE, name = "mem")
  if (any(membership[included] < 0 | membership[included] > 1, na.rm = TRUE)) {
    warning("Check if 'mem' values are degree of membership in a fuzzy classification.")
  }
  values <- rep(FALSE, terra::ncell(mem))
  groups <- split(which(included), labels$index[included])
  for (cells in groups) {
    cells <- cells[!is.na(membership[cells])]
    count <- max(0, min(length(cells), round(sum(membership[cells]))))
    if (count > 0) {
      ranked <- order(membership[cells], decreasing = TRUE, method = "radix")
      values[cells[ranked[seq_len(count)]]] <- TRUE
    }
  }
  bin <- terra::rast(mem)
  terra::values(bin) <- values
  bin
}
