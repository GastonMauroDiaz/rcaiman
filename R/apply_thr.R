#' Apply threshold
#'
#' Global or local thresholding of images.
#'
#' It is a wrapper function around the operator `>` from the `terra` package. If
#' a single threshold value is provided as the `thr` argument, it is applied to
#' every pixel of the object `r`. If instead a [SpatRaster-class] is provided,
#' a particular threshold is applied to each particular pixel.
#'
#' @param r [SpatRaster-class]. A greyscale image.
#' @param thr Numeric vector of length one or a single-layer raster from the
#'   class [SpatRaster-class]. Threshold.
#'
#' @return An object of class [SpatRaster-class] with values `0` and `1`.
#' @export
#'
#' @examples
#' r <- read_caim()
#' bin <- apply_thr(r$Blue, thr_isodata(r$Blue[]))
#' plot(bin)
#' \dontrun{
#' # This function is useful in combination with the ‘autothresholdr’
#' # package. For example:
#' require(autothresholdr)
#' thr <- auto_thresh(r$Blue[], "IsoData")[1]
#' bin <- apply_thr(r$Blue, thr)
#' plot(bin)
#' }
apply_thr <- function (r, thr)
{
  .is_single_layer_raster(r)
  r[is.na(r)] <- 0

  if (any(is(thr, "numeric"), is(thr, "integer"))) {
    stopifnot(length(thr) == 1)
    tmp <- values(r)
    if (thr < min(tmp, na.rm = TRUE))
      stop("\"thr\" should be greater than or equal to minimum layer value")
    if (thr >= max(tmp, na.rm = TRUE))
      stop("\"thr\" should be lower than maximum layer value")
  } else {
    if (!is(thr, "SpatRaster"))
      stop(paste("\"thr\" class should be \"numeric\",",
                 "\"integer\", or \"SpatRaster\""))
  }
  r > thr
}
