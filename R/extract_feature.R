#' Extract feature
#'
#' Extract features from raster images.
#'
#' Given a single-layer raster, a segmentation, and a function,
#' `extract_features` will return a numeric vector or a [SpatRaster-class]
#' depending on whether the parameter `return_raster` is `TRUE` or `FALSE`. For
#' the first case, each pixel of each segment will adopt the respective
#' extracted feature value. For the second case, the return will be a vector of
#' length equal to the total number of segments. Each value will be obtained by
#' processing all pixels that belong to a segment with the provided function.
#'
#' @param r [SpatRaster-class]. Single layer raster.
#' @param segmentation [SpatRaster-class]. The segmentation of `r`.
#' @param fun A `function` that takes a vector as input and returns a
#'   one-length numeric or logical vector as output (e.g. mean).
#' @param return_raster Logical vector of length one, see details.
#' @param ignore_label_0 Logical vector of length one. If this is `TRUE`,
#'   the segment labeled with `0` will be ignored.
#'
#' @export
#'
#' @return If `return_raster` is set to `TRUE`, then an object of
#'   class [SpatRaster-class] with the same pixel dimensions than `r`
#'   will be returned. Otherwise, the return is a numeric vector of length equal
#'   to the number of segments found in `segmentation`.
#'
#' @examples
#' r <- read_caim()
#' z <- zenith_image(ncol(r),lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 10)
#' print(extract_feature(r$Blue, g, return_raster = FALSE))
#' # plot(extract_feature(r$Blue, g, return_raster = TRUE))
extract_feature <- function(r, segmentation,
                            fun = mean,
                            return_raster = TRUE,
                            ignore_label_0 = TRUE) {

  .is_single_layer_raster(r)
  .is_single_layer_raster(segmentation, "segmentation")
  stopifnot(any(class(fun) == "function", class(fun) == "standardGeneric"))
  stopifnot(class(return_raster) == "logical")
  stopifnot(compareGeom(r, segmentation) == TRUE)

  if (ignore_label_0 == TRUE) segmentation[segmentation == 0] <- NA

  feature <- tapply(terra::values(r), terra::values(segmentation), fun)

  if (return_raster) {
    id <- as.numeric(names(feature))
    return(terra::subst(segmentation, id, feature))
  } else {
    ids <- names(feature)
    feature <- as.numeric(feature)
    names(feature) <- ids
    return(feature)
  }
}
