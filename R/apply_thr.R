#' Apply threshold
#'
#' Global or local thresholding of images.
#'
#' It is a wrapper function around the operator \code{>} from the ‘terra’
#' package. If a single threshold value is provided as \code{thr} argument, it
#' is applied to every pixel of the object \code{r}. If, instead, a
#' \linkS4class{SpatRaster} is provided as \code{thr} argument, then a
#' particular threshold is applied to each particular pixel.
#'
#' @param r \linkS4class{SpatRaster}. A greyscale image.
#' @param thr Numeric vector of length one or \linkS4class{SpatRaster}.
#'   Threshold.
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#' @export
#'
#' @family Binarization Functions
#'
#' @examples
#' r <- read_caim()
#' apply_thr(r$Blue, 120)
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
