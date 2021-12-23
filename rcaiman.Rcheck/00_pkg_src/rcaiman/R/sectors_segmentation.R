#' Sectors segmentation
#'
#' Segmenting an hemispherical view by slicing the azimuth angle from \code{0}
#' to \code{360º} in equals intervals.
#'
#' @param a \code{\linkS4class{RasterLayer}} built with
#'   \code{\link{azimuth_image}}.
#' @inheritParams rings_segmentation
#'
#' @return \linkS4class{RasterLayer} with segments shaped like pizza slices.
#' @export
#'
#' @family Segmentation functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' a <- azimuth_image(z)
#' sectors <- sectors_segmentation(a, 15)
#' plot(sectors == 1)
sectors_segmentation <- function(a, angle_width, return_angle = FALSE) {
  stopifnot(class(a) == "RasterLayer")
  stopifnot(class(return_angle) == "logical")
  stopifnot(length(angle_width) == 1)

  if (!.is_whole(360 / angle_width)) {
    stop(
      paste("angle_width should divide the 0 to 360",
            "range into a whole number of segments.")
    )
  }

  intervals <- seq(0, 360, angle_width)
  c1 <- intervals[1:(length(intervals) - 1)]
  c2 <- intervals[2:length(intervals)]
  if (return_angle) {
    c3 <- (c1 + c2) / 2
  } else {
    c3 <- 1:(length(intervals) - 1)
  }

  rcl <- matrix(c(c1, c2, c3), ncol = 3)
  reclassify(a, rcl)
}
