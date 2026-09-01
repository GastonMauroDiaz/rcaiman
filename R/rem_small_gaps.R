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
#'   area (\eqn{A}) and radii standard deviation (\eqn{\text{SD}_r}) of canopy gaps. While \eqn{A}
#'   is simply the total number of pixels of the segment, \eqn{\text{SD}_r} is the standard
#'   deviation of the Euclidean distance between the gap center and the pixels
#'   that are in the gap contour. The center is computed as the average of the
#'   coordinates of the contour pixels. The effective circular area equation is
#'   \eqn{\text{ECA} = A / SD_r}, providing a size metric that accounts for deviations from
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
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' thr <- apply_by_direction(r, z, a, m)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
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
