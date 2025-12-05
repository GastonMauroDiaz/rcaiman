#' Remove small gaps
#'
#' Remove patches of connected `TRUE` pixels based on a white top-hat pixel-based
#' morphological filter or an object-based fixed area/size threshold.
#'
#' @details
#' Methods provided:
#' \describe{
#'   \item{pixel-based:}{Internally runs [EBImage::whiteTopHat()] with a kernel
#'   of size equal to `size` and diamond shape.}
#'   \item{object-based:}{
#'   Patches of connected `TRUE` pixes (canopy gaps) are turn `FALSE` if they are
#'   rounded but too small, or if they are of a shape and size unlikely for a
#'   real canopy gap. The function internally uses _EBImage_ for calculates the
#'   area (\eqn{A}) and radii standard deviation (\eqn{SD_r}) of canopy gaps. While \eqn{A}
#'   is simply the total number of pixels of the segment, \eqn{SD_r} is the standard
#'   deviation of the Euclidean distance between the gap center and the pixels
#'   that are in the gap contour. The center is computed as the average of the
#'   coordinates of the contour pixels. The effective circular area equation is
#'   \eqn{ECA = A / SD_r}, providing a size metric that accounts for deviations from
#'   circularity. With this method, argument `size` is a threshold of ECA for
#'   recognizing small gaps.
#'   }
#' }
#'
#' @param size numeric vector of length one. See *Details*.
#' @param method character vector of length one. See *Details*.
#'
#' @inheritParams compute_canopy_openness
#'
#' @return Logical [terra::SpatRaster-class] with the same dimensions as `bin`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' r <- read_caim()
#' bin <- binarize_with_thr(r$Blue, thr_isodata(r$Blue[]))
#' plot(bin)
#' bin <- rem_small_gaps(bin, 5, method = "pixel-based")
#' plot(bin)
#' bin <- rem_small_gaps(bin, 11, method = "object-based")
#' plot(bin)
#' }
rem_small_gaps <- function(bin, size = 9, method = "pixel-based") {
  .this_requires_EBImage()
  .assert_logical_mask(bin)
  .check_vector(size, "numeric", 1, sign = "positive")
  .assert_choice(method, c("pixel-based", "object-based"))

  if (method == "pixel-based") {
    kern <-  EBImage::makeBrush(size, shape = "diamond")
    small_gaps <-  terra::setValues(bin,
                                    EBImage::whiteTopHat(as.array(bin), kern))
    bin[small_gaps] <- 0
  } else {
    bwlabels <- EBImage::bwlabel(as.array(bin))
    shape <- EBImage::computeFeatures.shape(bwlabels)
    bwlabels <- terra::setValues(bin, bwlabels)
    shape <- terra::subst(bwlabels, 1:nrow(shape), shape)

    ## Use effective circular area (ECA)
    eca <- shape$s.area / (shape$s.radius.sd + 1)
    bin[eca < size] <- 0
  }
  as.logical(bin)
}
