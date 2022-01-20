#' Rings segmentation
#'
#' Segmenting an hemispherical view by slicing the zenith angle from \code{0} to
#' \code{90º} in equals intervals.
#'
#' @inheritParams azimuth_image
#' @param angle_width Numeric vector of length one. Angle in degrees able to
#'   divide the angle range into a whole number of segments.
#' @param return_angle Logical vector of length one. If it is \code{FALSE}, all
#'   the pixels that belong to a segment are labeled with an ID number.
#'   Otherwise, the angle mean of the segment is assigned to the pixels.
#'
#' @return An object from the class \linkS4class{RasterLayer} with segments
#'   shaped like concentric rings.
#' @export
#'
#' @family Segmentation functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' rings <- rings_segmentation(z, 15)
#' plot(rings == 1)
rings_segmentation <- function(z, angle_width, return_angle = FALSE) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(return_angle) == "logical")
  stopifnot(length(angle_width) == 1)

  if (!.is_whole(90 / angle_width)) {
    stop(
      paste("angle_width should divide the,",
            "0 to 90 range into a whole number of segments.")
    )
  }

  intervals <- seq(0, 90, angle_width)
  c1 <- intervals[1:(length(intervals) - 1)]
  c2 <- intervals[2:length(intervals)]

  if (return_angle) {
    c3 <- (c1 + c2) / 2
  } else {
    c3 <- 1:(length(intervals) - 1)
  }
  rcl <- matrix(c(c1, c2, c3), ncol = 3)
  reclassify(z, rcl)
}
