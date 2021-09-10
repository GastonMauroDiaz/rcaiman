#' Apply a global threshold
#'
#' Apply a global threshold value to an image.
#'
#' It is a wrapper function around this operator from the 'raster' package:
#' \code{>}. The same threshold value is applied to every layer of the raster
#' object. The most frequent use is ingesting a single layered object, often
#' called a greyscale images.
#'
#' @param r \linkS4class{Raster}.
#' @param thr Numeric vector of length one. Threshold value.
#'
#' @return \linkS4class{RasterLayer} with values \code{0} and \code{1}.
#' @export
#'
#' @examples
#' r <- read_caim()
#' apply_thr(r$Blue, 120)
#' \dontrun{
#' # This function is very useful in combination with the 'autothresholdr'
#' package. For examples:
#' thr <- autothresholdr::auto_thresh(r$Blue[], "IsoData")
#' apply_thr(r$Blue, thr)
#' }
apply_thr <- function (r, thr)
{
  stopifnot(any(class(thr) == "numeric", class(thr) == "integer"))
  stopifnot(length(thr) == 1)

  tmp <- values(r)
  if (thr < min(tmp, na.rm = TRUE))
    stop("thr must be greater than or equal to minimum layer value")
  if (thr >= max(tmp, na.rm = TRUE))
    stop("thr must be lower than maximum layer value")

  r > thr

}
