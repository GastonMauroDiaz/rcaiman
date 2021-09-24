#' Apply threshold
#'
#' Global or local thresholding of images.
#'
#' It is a wrapper function around this operator from the 'raster' package:
#' \code{>}.
#'
#' If a single threshold value is provided as \code{thr} argument, it is applied
#' to every pixel of the raster object \code{r}.
#'
#' If a \linkS4class{RasterLayer} is provided as \code{thr}, then the operation
#' is performed pixel by pixel, so it allows local threholding. In other words,
#' each pixel can have a personalized threshold.
#'
#' @param r \linkS4class{RasterLayer}.
#' @param thr Numeric vector of length one or \linkS4class{RasterLayer}.
#'   Threshold value or values.
#'
#' @return \linkS4class{RasterLayer} with values \code{0} and \code{1}.
#' @export
#'
#' @seealso \code{\link{regional_thresholding}}.
#'
#' @examples
#' r <- read_caim()
#' apply_thr(r$Blue, 120)
#' \dontrun{
#' # This function is very useful in combination with the 'autothresholdr'
#' package. For examples:
#' require(autothresholdr)
#' thr <- auto_thresh(r$Blue[], "IsoData")[1]
#' apply_thr(r$Blue, thr)
#' }
apply_thr <- function (r, thr)
{

  stopifnot(class(r) == "RasterLayer" )

  if (any(class(thr) == "numeric", class(thr) == "integer")) {
    stopifnot(length(thr) == 1)
    tmp <- values(r)
    if (thr < min(tmp, na.rm = TRUE))
      stop("\"thr\" should be greater than or equal to minimum layer value")
    if (thr >= max(tmp, na.rm = TRUE))
      stop("\"thr\" should be lower than maximum layer value")
  } else {
    if (class(thr) != "RasterLayer")
      stop(paste("\"thr\" class should be \"numeric\",",
                 "\"integer\", or \"RasterLayer\""))
  }

  r > thr
}
