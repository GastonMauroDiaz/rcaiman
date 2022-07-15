#' Sectors segmentation
#'
#' Segmenting an hemispherical view by slicing the azimuth angle from zero
#' to 360º in equals intervals.
#'
#' @inheritParams ootb_mblt
#' @inheritParams rings_segmentation
#'
#' @return An object from the class \linkS4class{SpatRaster} with segments
#'   shaped like pizza slices.
#' @export
#'
#' @family Segmentation Functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' a <- azimuth_image(z)
#' sectors <- sectors_segmentation(a, 15)
#' plot(sectors == 1)
sectors_segmentation <- function(a, angle_width, return_angle = FALSE) {
  .is_single_layer_raster(a, "a")
  stopifnot(.get_max(a) <= 360)
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
  sectors <- terra::classify(a, rcl)
  sectors[is.na(sectors)] <- 0
  names(sectors) <- "Sectors"
  sectors
}
