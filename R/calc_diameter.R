#' Calculate diameter
#'
#' Calculate the diameter in pixels of a 180º fisheye image.
#'
#' This function helps handle devices with a field of view different than 180º.
#' Given a lens projection function and data points consisting of radii
#' (pixels) and their correspondent zenith angle (\eqn{\theta}), it returns the
#' horizon radius (i.e., the radius for \eqn{\theta} equal to 90º).
#'
#' When working with non-circular hemispherical photography, this function will
#' help to find the diameter that a circular image would have if the equipment
#' would record the whole hemisphere.
#'
#' The required data (radius-angle data pairs) can be obtained following the
#' instructions given in the
#' [user manual of Hemisfer software](https://www.schleppi.ch/patrick/hemisfer/help/en/lens.htm).
#' The following is a slightly simpler alternative:
#'
#' * Find a vertical wall and a leveled floor, both well-constructed.
#'
#' * Draw a triangle of \eqn{5 \times 4 \times 3} meters on the floor, with
#' the 4-meter side over the wall.
#'
#' * Locate the camera over the vertex that is 3 meters away from the wall.
#' Place it at a given height above the floor, 1.3 meters for instance.
#'
#' * Make a mark on the wall at the chosen height over the wall-vertex nearest
#' to the camera-vertex. Make four more marks with one meter of spacing and
#' following a horizontal line. This will create marks for 0º, 18º, 34º, 45º,
#' and 54º \eqn{\theta}.
#'
#' * Before taking the photograph, do not forget to align the zenith
#' coordinates with the 0º \eqn{\theta} mark and check if the optical axis is
#' leveled.
#'
#' The [line selection
#' tool](https://imagej.nih.gov/ij/docs/guide/146-19.html#toc-Subsection-19.2)
#' of [ImageJ](https://imagej.nih.gov/ij/download.html) can be used to measure
#' the distance in pixels between points on the image. Draw a line, and use the
#' dropdown menu Analyze>Measure to obtain its length.
#'
#' For obtaining the projection of a new lens, refer to
#' [calibrate_lens()].
#'
#'
#' @inheritParams zenith_image
#' @param radius Numeric vector. Distance in pixels from the zenith.
#' @param angle Numeric vector. Zenith angle in degrees.
#'
#' @family Lens Functions
#'
#' @return Numeric vector of length one. Diameter adjusted to a whole number
#'   (see [zenith_image()] for details about that constrain).
#'
#' @export
#'
#' @examples
#' # Nikon D50 and Fisheye Nikkor 10.5mm lens
#' calc_diameter(lens("Nikkor_10.5mm"), 1202, 54)
calc_diameter <- function(lens_coef, radius, angle) {
  stopifnot(length(radius) == length(angle))

  Rfor90 <- calc_relative_radius(90, lens_coef)
  RforMyAngle <- calc_relative_radius(angle, lens_coef)

  fun <- function(radius, RforMyAngle) {
    Rfor90 * radius / RforMyAngle * 2
  }

  if (length(radius) == 1) {
    diameter <- round(fun(radius, RforMyAngle))
  } else {
    diameters <- unlist(Map(fun, radius, RforMyAngle))
    diameter <- round(stats::median(diameters))
    attr(diameter, "IQR") <- stats::IQR(diameters)
  }

  if (!.is_even(diameter)) {
    diameter <- diameter + 1
  }

  diameter
}

