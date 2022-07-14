#' Calculate zenith raster coordinates
#'
#' Calculate zenith raster coordinates from points digitized with the
#' open-source software package ‘ImageJ’. The zenith is the point on the image
#' that represents the zenith when upward-looking photographs are taken with the
#' optical axis parallel to the vertical line.
#'
#' The technique described under the headline ‘Optical center characterization’
#' of the
#' \href{https://www6.paca.inrae.fr/can-eye/content/download/3052/30819/version/4/file/CAN_EYE_User_Manual.pdf}{user
#' manual of the software Can-Eye} can be used to acquire the data for
#' determining the zenith coordinates. This technique was used by
#' \insertCite{Pekin2009;textual}{rcaiman}, among others. Briefly, it consists
#' in drilling a small hole in the cap of the fisheye lens (it must be away from
#' the center of the cap), and taking about ten photographs without removing the
#' cap. The cap must be rotated about 30º before taking each photograph.
#' \strong{The method implemented here do not support multiple holes}.
#'
#' The
#' \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#sec:Multi-point-Tool}{point
#' selection tool of ‘ImageJ’ software} should be used to manually digitize the
#' white dots and create a CSV file to feed this function.
#'
#' Another method --only valid for circular hemispherical photographs-- is
#' taking a very bright picture (for example, a picture of a room with walls
#' painted in light colors) with the lens completely free (do not use any
#' mount). Then, digitize points over the perimeter of the circle. This was the
#' method used for producing the example file (see Examples). It is worth noting
#' that the perimeter of the circle depicted in a circular hemispherical
#' photograph is not necessarily the horizon.
#'
#' @inheritParams calibrate_lens
#'
#' @references \insertAllCited{}
#'
#' @family Lens Functions
#'
#' @export
#'
#' @return Numeric vector of length two. Raster coordinates of the zenith,
#'   assuming a lens facing up with its optical axis parallel to the vertical
#'   line. It is important to note the difference between the raster coordinates
#'   and the Cartesian coordinates. In the latter,  the vertical axis value
#'   decreases downward, but the opposite is true for the raster coordinates,
#'   which works like a spreadsheet.
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/points_over_perimeter.csv",
#'                     package = "rcaiman")
#' calc_zenith_raster_coord(path)
#' }
calc_zenith_raster_coord <- function(path_to_csv) {
  if (!requireNamespace("conicfit", quietly = TRUE)) {
    stop(paste("Package \"conicfit\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
  requireNamespace("conicfit", quietly = TRUE)
  x <- utils::read.csv(path_to_csv)[, -(1:5)]
  circle <- as.matrix(x) %>%
      conicfit::CircleFitByKasa() %>%
      .[-3]
  names(circle) <- c("col", "row")
  circle
}
