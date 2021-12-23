#' Apply threshold
#'
#' Global or local thresholding of images.
#'
#' It is a wrapper function around the operator \code{>} from the ‘raster’
#' package. If a single threshold value is provided as \code{thr} argument, it
#' is applied to every pixel of the raster object \code{r}. If instead a
#' \linkS4class{RasterLayer} is provided, then a particular threshold
#' is applied to each particular pixel.
#'
#' @param r \linkS4class{RasterLayer}
#' @param thr Numeric vector of length one or \linkS4class{RasterLayer}.
#'   Threshold.
#'
#' @return \linkS4class{RasterLayer} with values \code{0} and \code{1}.
#' @export
#'
#' @family Tools functions
#'
#' @seealso \code{\link{regional_thresholding}}.
#'
#' @examples
#' r <- read_caim()
#' apply_thr(r$Blue, 120)
#' \dontrun{
#' # This function is useful in combination with the ‘autothresholdr’
#' package. For examples:
#' require(autothresholdr)
#' thr <- auto_thresh(r$Blue[], "IsoData")[1]
#' bin <- apply_thr(r$Blue, thr)
#' plot(bin)
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
