#' Binarize with known thresholds
#'
#' @description
#' Apply a threshold values or a raster with threshold values to a grayscale
#' image, producing a binary image.
#'
#' @details
#' This function supports both global and pixel-wise thresholding. It is a
#' wrapper around the `>` operator from the `terra` package. If a single numeric
#' threshold is provided via `thr`, it is applied globally to all pixels in `r`.
#' If instead a [terra::SpatRaster-class] object is provided, local thresholding
#' is performed, where each pixel is compared to its corresponding threshold
#' value.
#'
#' This is useful after estimating thresholds using [thr_two_corner()],
#' [thr_isodata()], or
#' [apply_by_direction()] with `fun = thr_isodata`, among other possibilities.
#'
#' @param r numeric [terra::SpatRaster-class] of one layer.
#' @param thr either a numeric vector of length one (for global thresholding) or
#'   a numeric [terra::SpatRaster-class] with one layer (for local thresholding).
#'
#' @return Logical [terra::SpatRaster-class] (`TRUE` for sky, `FALSE` for
#'   non-sky) with the same dimensions as `r`.
#'
#' @note
#' For global thresholding, `thr` must be greater than or equal to the minimum
#' value of `r` and lower than its maximum value.
#'
#' @export
#'
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, fill_na = "none")$Blue
#' # Use a small crop to illustrate thresholding.
#' r <- crop_caim(r, c(600, 800), 200, 200)
#' bin <- binarize_with_thr(r, thr_isodata(r[]))
#' plot(bin)
#'
#' \dontrun{
#' # This function is also compatible with thresholds estimated using
#' # the 'autothresholdr' package:
#' require(autothresholdr)
#' r <- normalize_min_max(r) %>% multiply_by(255) %>% round()
#' thr <- auto_thresh(r[], "IsoData")[1]
#' bin <- binarize_with_thr(r, thr)
#' plot(bin)
#' }
binarize_with_thr <- function (r, thr)
{
  .assert_single_layer(r)

  if (is.list(thr)) {
    stop(sprintf("`%s` cannot be a list.", "thr"),
         call. = FALSE)
  }

  if (tryCatch(
    .check_vector(thr, "numeric", sign = "any"),
    error = function(e)
      FALSE
  )) {
    if (length(thr) != 1) {
      stop("`thr` must be of length one.")
    }
    v <- values(r)
    if (thr < min(v, na.rm = TRUE))
      stop("`thr` must be greater than or equal to minimum layer value.")
    if (thr >= max(v, na.rm = TRUE))
      stop("`thr` must be lower than maximum layer value.")
  } else {
    .assert_single_layer(thr)
  }
  bin <- r > thr
  bin[is.na(bin)] <- FALSE
  bin
}
