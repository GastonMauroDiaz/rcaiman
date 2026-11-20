#' Assign zenith-ring labels
#'
#' @description
#' Divide the sky into rings of equal zenith-angle width and label the
#' observed directions on the image grid.
#'
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, in that order. Only `zenith` is used;
#'   `azimuth` is ignored, including its missing values.
#' @param angle_width positive numeric vector of length one. Ring width in
#'   degrees. Must divide 90 into a whole number of rings.
#' @param return character vector of length one. Output mode: `"id"` (default)
#'   or `"angle"`.
#'
#' @details
#' Rings are defined relative to the upward vertical, independently of the
#' optical axis. Supply the acquisition orientation to [direction_image()]
#' before segmenting. A partial view retains the same ring boundaries and IDs
#' as a complete view of the sky.
#'
#' Intervals include their upper boundary; the first also includes zero.
#' For example, with a width of 15 deg, the first ring includes 0 through
#' 15 deg, and the second includes angles greater than 15 through 30 deg.
#' Missing zenith angles and directions below the horizon receive zero.
#' A direction exactly at zenith belongs to the first ring even when its
#' azimuth is undefined.
#'
#' @return Single-layer [terra::SpatRaster-class] on the same grid as `d`.
#'   Values are ring IDs starting at 1, or angular midpoints in degrees when
#'   `return = "angle"`. Zero denotes cells outside the segmentation.
#'   The object carries attribute `return`.
#' @seealso [sector_segmentation()], [sky_grid_segmentation()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' rings <- ring_segmentation(d, 15)
#' plot(rings)
ring_segmentation <- function(d, angle_width, return = "id") {
  z <- .sky_segmentation_layer(d, "zenith")
  n <- .sky_interval_count(angle_width, 90)
  .assert_choice(return, c("id", "angle"))
  bounds <- seq(0, 90, length.out = n + 1)
  labels <- if (return == "angle") {
    (head(bounds, -1) + tail(bounds, -1)) / 2
  } else seq_len(n)
  rings <- terra::classify(z, cbind(head(bounds, -1), tail(bounds, -1), labels),
                           right = TRUE, include.lowest = TRUE, others = 0)
  rings[is.na(rings)] <- 0
  names(rings) <- "Rings"
  attr(rings, "return") <- return
  rings
}

.sky_segmentation_layer <- function(d, layer) {
  .assert_spatraster(d)
  if (!identical(names(d), c("zenith", "azimuth"))) {
    stop("`d` must have two layers named zenith and azimuth, in that order.",
         call. = FALSE)
  }
  angle <- terra::values(d[[layer]], mat = FALSE)
  .check_vector(angle, "numeric", allow_na = TRUE, name = paste0("d$", layer))
  invalid <- if (layer == "zenith") angle < 0 | angle > 180 else
    angle < 0 | angle >= 360
  if (any(invalid, na.rm = TRUE)) {
    stop("`d$", layer, "` must be in ",
         if (layer == "zenith") "[0, 180]" else "[0, 360)",
         " degrees.", call. = FALSE)
  }
  d[[layer]]
}

.sky_interval_count <- function(angle_width, extent) {
  .check_vector(angle_width, "numeric", 1, sign = "positive")
  n <- extent / angle_width
  if (!is.finite(n) || n < 1 ||
      abs(n - round(n)) > 64 * .Machine$double.eps * max(1, n)) {
    stop("`angle_width` must divide ", extent,
         " into a whole number of intervals.", call. = FALSE)
  }
  round(n)
}
