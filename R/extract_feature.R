#' Extract feature
#'
#' Extract features from raster images.
#'
#' Given a single layer raster, a segmentation, and a function,
#' \code{extract_features} will returns a numeric vector or a
#' \linkS4class{RasterLayer} depending on whether the parameter
#' \code{return_raster} is \code{TRUE} or \code{FALSE}. For the first case, each
#' pixel of each segment will adopt the respective extracted feature value. For
#' the second case, the return will be the extracted feature as a vector of
#' length equal to the total number of segments. Each extracted feature value
#' will be obtained by processing all pixels that belong to a segment with the
#' provided function.
#'
#' @param r \linkS4class{RasterLayer}. Single layer raster.
#' @param segmentation \linkS4class{RasterLayer}. The segmentation of \code{r}.
#' @param fun \code{function} that takes a vector as input and returns a
#'   one-length numeric or logical vector as output (e.g. mean).
#' @param return_raster Logical vector of length one.
#'
#' @family Tools functions
#'
#' @export
#'
#' @return If \code{return_raster} is set to \code{TRUE}, then an object of
#'   class \linkS4class{RasterLayer} with the same pixel dimensions than
#'   \code{r} will be returned. Otherwise, the return is a numeric vector of
#'   length equal to the number of segments found in \code{segmentation}.
#'
#' @examples
#' \dontrun{
#' r <- read_caim()
#' z <- zenith_image(ncol(r),lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' print(extract_feature(r$Blue, g, return_raster = FALSE))
#' plot(extract_feature(r$Blue, g, return_raster = TRUE))
#' }
extract_feature <- function(r, segmentation, fun = mean, return_raster = TRUE) {

  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(segmentation) == "RasterLayer")
  stopifnot(any(class(fun) == "function", class(fun) == "standardGeneric"))
  stopifnot(class(return_raster) == "logical")

  feature <- tapply(values(r), values(segmentation), fun)

  if (return_raster) {
    id <- as.numeric(names(feature))
    df <- data.frame(id, feature)
    subs(segmentation, df)
  } else {
    ids <- names(feature)
    feature <- as.numeric(feature)
    names(feature) <- ids
    feature
  }
}
