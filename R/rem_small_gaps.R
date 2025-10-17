#' Remove small gaps
#'
#' Remove patches of connected `TRUE` pixels based on fixed area/size
#' thresholds.
#'
#' Pathes of connected `TRUE` pixes (canopy gaps) are turn `FALSE` if they are
#' rounded but too small, or if they are of a shape and size unlikely for a real
#' canopy gap. The function internally uses _EBImage_ for calculates the area
#' ($A$) and radii standard deviation ($SD_r$) of canopy gaps. While $A$ is
#' simply the total number of pixels of the segment, $SD_r$ is the standard
#' deviation of the Euclidean distance between the gap center and the pixels
#' that are in the gap contour. The center is computed as the average of the
#' coordinates of the contour pixels. The effective circular area equation is
#' $ECA = A / SD_r$, providing a size metric that accounts for deviations from
#' circularity.
#'
#' @param eca_thr numeric vector of length one. Threshold of effective circular
#'   area for recognizing small gaps.
#'
#' @inheritParams compute_canopy_openness
#'
#' @return Logical [terra::SpatRaster-class] with the same dimensions as `bin`.
#' Compared to the input `bin`, black regions (`FALSE`) have been expanded by
#' the specified buffer distance.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' r <- read_caim()
#' bin <- binarize_with_thr(r$Blue, thr_isodata(r$Blue[]))
#' plot(bin)
#' bin <- grow_black(bin, 11)
#' plot(bin)
#' }
rem_small_gaps <- function(bin, eca_thr = 9) {
  .this_requires_EBImage()
  .assert_logical_mask(bin)
  .check_vector(eca_thr, "numeric", 1, sign = "positive")

  bwlabels <- EBImage::bwlabel(as.array(bin))
  shape <- EBImage::computeFeatures.shape(bwlabels)
  bwlabels <- terra::setValues(bin, bwlabels)
  shape <- terra::subst(bwlabels, 1:nrow(shape), shape)

  ## Use effective circular area (ECA)
  eca <- shape$s.area / (shape$s.radius.sd + 1)
  bin[eca < eca_thr] <- 0
  bin
}
