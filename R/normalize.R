normalize <- function(r, mn, mx, ...) {
  stopifnot(length(mn) == 1)
  stopifnot(length(mx) == 1)
  (r - mn) / (mx - mn)
}

#' Normalize data
#'
#' Normalize data lying between \code{mn} and \code{mx} in the range \code{0} to
#' \code{1}. Data greater than \code{mx} get values greater than \code{1} in a
#' proportional fashion. Conversely, data less than \code{mn} get values less
#' than \code{0}.
#'
#' @param r \code{\linkS4class{Raster}}.
#' @param mn Numeric vector of length one. Minimum expected value.
#' @param mx Numeric vector of length one. Maximum expected value.
#' @param ... Additional arguments as for \code{\link[raster]{writeRaster}}.
#'
#' @examples
#' normalize(read_caim(), 0, 255)
#'
setGeneric("normalize", normalize)
#' @export normalize

#' @describeIn normalize Result is a \linkS4class{Raster} of the same type that
#'   \code{r}.
setMethod("normalize",
          signature(r = "Raster"),
          function(r, mn, mx, ...) {
            fun <- .makeF8multi(function(r,...) normalize(r, mn, mx), ...)
            fun(r, ...)
          }
)

#' @rdname normalize
setMethod("normalize",
          signature(r = "RasterLayer"),
          function(r, mn, mx, ...) {
            fun <- .makeF8single(function(r,...) normalize(r, mn, mx), ...)
            fun(r, ...)
          }
)
