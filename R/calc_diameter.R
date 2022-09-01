#' Calculate diameter
#'
#' Calculate the diameter in pixels of a 180º fisheye image.
#'
#' This function helps handle devices with a field of view different than 180
#' degrees. Given a lens projection function and data points consisting of radii
#' (pixels) and their correspondent zenith angle (\eqn{\theta}), it returns the
#' radius of the horizon (i.e., the radius for \eqn{\theta} equal to 90
#' degrees).
#'
#' When working with non-circular hemispherical photography, this function will
#' help to find the diameter that a circular image would have if the equipment
#' would record the whole hemisphere.
#'
#' The required data (radius-angle data) can be obtained following the
#' instructions given in the
#' \href{https://www.schleppi.ch/patrick/hemisfer/help/en/lens.htm}{user manual
#' of Hemisfer software}. The following is a slightly simpler alternative:
#'
#' \itemize{
#'
#' \item Find a vertical wall and a leveled floor, both well-constructed. For
#' instance, a parking lot.
#'
#' \item Draw a triangle of \eqn{5 \times 4 \times 3} meters on the floor, with
#' the 4-meter side over the wall.
#'
#' \item Locate the camera over the vertex that is 3 meters away from the wall.
#' Place it at a given height above the floor, 1.3 meters for instance.
#'
#' \item Make a mark on the wall at 1.3 meters over the wall-vertex nearest to
#' the camera-vertex. Make four more marks with one meter of spacing and
#' following a horizontal line. This will create marks for 0º, 18º, 34º, 45º,
#' and 54º \eqn{\theta}.
#'
#' \item Before taking the photograph, don not forget to align the zenith
#' coordinates with the 0º \eqn{\theta} mark and check if the optical axis is
#' leveled.
#'
#' }
#'
#' The
#' \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#toc-Subsection-19.2}{line
#' selection tool} of \href{https://imagej.nih.gov/ij/download.html}{ImageJ} can
#' be used to measure the distance in pixels between points on the image. Draw a
#' line, and use the dropdown menu Analyze>Measure to obtain its length.
#'
#' For obtaining the projection of a new lens, refer to
#' \code{\link{calibrate_lens}}.
#'
#'
#' @inheritParams zenith_image
#' @param radius_px Numeric vector. Distance in pixels from the zenith.
#' @param angle Numeric vector. Zenith angle in degrees.
#'
#' @family Lens Functions
#'
#' @return Numeric vector of length one. Diameter adjusted to a whole number
#'   (see \code{\link{zenith_image}} for details about that constrain).
#'
#' @export
#'
#' @examples # Nikon D50 and Fisheye Nikkor 10.5 mm lens
#' calc_diameter(lens("Nikkor_10.5_mm"), 1202, 54)
calc_diameter <- function(lens_coef, radius_px, angle) {
  stopifnot(length(radius_px) == length(angle))

  Rfor90 <- calc_relative_radius(90, lens_coef)
  RforMyAngle <- calc_relative_radius(angle, lens_coef)

  fun <- function(radius_px, RforMyAngle) {
    Rfor90 * radius_px / RforMyAngle * 2
  }

  if (length(radius_px) == 1) {
    diameter <- round(fun(radius_px, RforMyAngle))
  } else {
    diameters <- unlist(Map(fun, radius_px, RforMyAngle))
    diameter <- round(stats::median(diameters))
    attr(diameter, "IQR") <- stats::IQR(diameters)
  }

  if (!.is_even(diameter)) {
    diameter <- diameter + 1
  }

  diameter
}

