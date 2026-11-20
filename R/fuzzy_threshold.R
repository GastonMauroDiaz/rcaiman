#' Fuzzy membership above a threshold
#'
#' Calculate the degree to which each pixel exceeds a supplied threshold.
#' The threshold and the width of the transition can vary independently
#' across the image.
#'
#' @param r Numeric single-layer [terra::SpatRaster-class].
#' @param thr Threshold in the units of `r`: a finite number or a single-layer
#'   raster on the grid of `r`, for example a local threshold image.
#' @param scale Logistic scale in the units of `r`: a finite nonnegative
#'   number or a single-layer raster on the grid of `r`. Larger values give
#'   a more gradual transition. Zero requests a sharp threshold.
#'
#' @details For positive scales, membership is
#'   `plogis((r - thr) / scale)`. It is 0.5 at the threshold, and rises
#'   from 0.25 to 0.75 between `thr - log(3) * scale` and
#'   `thr + log(3) * scale`. With zero scale, membership is 1 where
#'   the pixel exceeds the threshold and 0 elsewhere, including ties.
#'
#'   A spatially varying scale was used by
#'   \insertCite{Diaz2015;textual}{rcaiman} to adjust the sharpness of a
#'   brightness threshold using color information.
#'
#' @return A numeric raster on the grid of `r`, named `membership`, with
#'   values in `[0, 1]`. A missing value in any input raster gives a missing
#'   membership at that pixel.
#' @seealso [binarize_with_thr()], [thr_isodata()], [defuzzify()]
#' @references \insertAllCited{}
#' @export
#' @examples
#' r <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 5))
#' terra::values(r) <- rep(c(20, 30, 40, 50, 60), 2)
#' thr <- thr_isodata(r[])
#' mem <- fuzzy_threshold(r, thr, scale = 5)
#'
#' # Move the threshold in the second row; keep the same transition width.
#' local_thr <- terra::rast(r)
#' terra::values(local_thr) <- rep(c(thr, thr + 10), each = 5)
#' fuzzy_threshold(r, local_thr, scale = 5)[]
#'
#' # Keep the threshold fixed; make the transition sharp in the first row.
#' local_scale <- terra::rast(r)
#' terra::values(local_scale) <- rep(c(0, 5), each = 5)
#' fuzzy_threshold(r, thr, scale = local_scale)[]
fuzzy_threshold <- function(r, thr, scale) {
  .assert_single_layer(r)
  x <- terra::values(r, mat = FALSE)
  .check_vector(x, "numeric", allow_na = TRUE, name = "r")

  parameter <- function(value, name, sign) {
    if (inherits(value, "SpatRaster")) {
      .assert_single_layer(value, name = name)
      .assert_same_geom(r, value, y_name = name)
      value <- terra::values(value, mat = FALSE)
      .check_vector(value, "numeric", sign = sign, allow_na = TRUE, name = name)
    } else {
      .check_vector(value, "numeric", 1, sign = sign, name = name)
    }
    value
  }
  threshold <- parameter(thr, "thr", "any")
  width <- parameter(scale, "scale", "nonnegative")
  delta <- x - threshold
  sharp <- !is.na(width) & width == 0
  width[sharp] <- 1  # Evaluate the sharp cells separately from the logistic.
  membership <- stats::plogis(delta / width)
  membership[sharp] <- as.numeric(delta[sharp] > 0)

  out <- terra::rast(r)
  terra::values(out) <- membership
  names(out) <- "membership"
  out
}
