#' Calibrate lens
#'
#' Calibrate a fisheye lens. This type of lens has a wide field of view and
#' consistent azimuthal distortion. The latter property allows fitting a precise
#' mathematical relationship between the distance to the zenith on the image
#' space and the zenith angle on the hemispherical space.
#'
#' These are the instructions to produce the CSV file required by this function.
#' The following materials are required:
#'
#' \itemize{
#'
#' \item this package and \href{https://imagej.nih.gov/ij/download.html}{ImageJ}
#'
#' \item camera and lens
#'
#' \item tripod
#'
#' \item standard yoga mat
#'
#' \item table about two times wider than the yoga mat width
#'
#' \item twenty two push pins of different colors
#'
#' \item scissors
#'
#' \item one print of this \href{https://osf.io/tudzc/download}{sheet} (A1 size,
#' almost like a poster).
#'
#' }
#'
#' Cut the sheet by the dashed line. Place the yoga mat extended on top of the
#' table. Place the sheet on top of the yoga mat. Align the dashed line with the
#' yoga mat border closest to you. Place push pins on each cross. If you are
#' gentle, the yoga mat will allow you to do that without damaging the table.
#' Of course, other materials could be used to obtain the same result, such as
#' cardboard, foam, nails, etc.
#'
#' Place the camera on the tripod. Align its optical axis with the table while
#' looking for getting an image showing the overlapping of the three pairs of
#' push pins, as instructed in the print. In order to take care of the line of
#' pins at 90º relative to the optical axis, it may be of help to use the naked
#' eye to align the front of the lens with the pins (Strictly speaking, we need
#' to alight the nodal point of the lens instead of its front. The term
#' "entrance pupil" is also used to refer to this point, but least-parallax
#' point may be the best term).
#'
#' Transfer the photograph to the computer, open it with ImageJ, and use the
#' \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#sec:Multi-point-Tool}{point
#' selection tool} to digitize the push pins, starting from the zenith push pin
#' and not skipping any shown push pin. Then, use the dropdown menu
#' Analyze>Measure to open the window Results. To obtain the CSV, use File>Save
#' As...
#'
#' This method was inspired by the calibration board from
#' \insertCite{Clark1988;textual}{rcaiman}.
#'
#' \strong{TIP:} use \code{\link{test_lens_coef}} to test if coefficients are
#' OK. If not, try moving the last points a little bit. Putting the one of the
#' last push pin a few pixels farther from the zenith is usually enough. An
#' alternative is to round the coefficients, or truncate the last number of the
#' last coefficient.
#'
#' Consult this \href{https://docs.google.com/document/d/178yZDAcfx--Xn1Ye8Js-kUXuPCuYOHQL5fxAH7KBEoY/edit?usp=sharing}{document}
#' for additional details.
#'
#' @param path_to_csv Character vector of length one. Path to a CSV file created
#'   with the
#'   \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#sec:Multi-point-Tool}{point
#'    selection tool of ‘ImageJ’ software}.
#' @param degree Numeric vector of length one. Polynomial model degree.
#'
#' @return An object of class \emph{list} with named elements. \emph{lens_coef}
#'   stands for lens coefficients, \emph{max_theta} for maximum zenith angle in
#'   degrees, and \emph{max_theta_px} for distance in pixels between the zenith
#'   and the maximum zenith angle in pixels units. The latter should be
#'   double-checked, particularly if the zenith push pin is not exactly on the
#'   zenith pixel. To that end, do the following on ImageJ: use the
#'   \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#toc-Subsection-19.1}{rectangular
#'    selection tool} to create a small rectangle, open the Specify window by
#'   going to the dropdown menu Edit>Selection>Specify..., insert the zenith
#'   coordinates (obtained with \code{\link{calc_zenith_raster_coord}}) into the
#'   respective X and Y fields in order to align the upper-left corner of the
#'   rectangle with the zenith, mark it with the
#'   \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#toc-Subsection-19.14}{brush},
#'    use the
#'   \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#toc-Subsection-19.2}{straight
#'   selection tool} to find the length within the zenith and the maximum zenith
#'   angle showed in the image.
#'
#' @family Lens Functions
#'
#' @references \insertAllCited{}
#'
#'
#' @export
#'
#' @examples
#' path <- system.file("external/Results_calibration.csv", package = "rcaiman")
#' calibration <- calibrate_lens(path)
#' calibration$lens_coef
#' calibration$max_theta
#' calibration$max_theta_px
#' test_lens_coef(calibration$lens_coef)
calibrate_lens <- function(path_to_csv, degree = 3) {

  csv <- utils::read.csv(path_to_csv)
  csv <- cbind(csv$X, csv$Y)

  # center in (0,0)
  csv[, 1] <- csv[, 1] - csv[1, 1]
  csv[, 2] <- csv[, 2] - csv[1, 2]

  csv <- pracma::cart2pol(csv)
  px <- csv[, 2]
  max_fov_px <- px[length(px)]
  # remove the data point that is only useful for max_fox_px calculation
  px <- px[-length(px)]

  theta <- .degree2radian(seq(0, 90, 5))
  theta <- theta[1:length(px)] # because some lens would "see" less pushpins

  fit <- lm(theta ~ poly(px, degree, raw = TRUE) - 1)
  lens_coef_caneye <- unname(coefficients(fit))
  max_fov <- stats::predict(fit, data.frame(px = max_fov_px)) %>%
    .radian2degree()

  fit <- lm(px ~ poly(theta, degree, raw = TRUE) - 1)
  px90 <- stats::predict(fit, data.frame(theta = pi / 2))

  R <- px / px90
  fit <- lm(R ~ poly(theta, degree, raw = TRUE) - 1)

  list(lens_coef = unname(coefficients(fit)),
       lens_coef_caneye = lens_coef_caneye,
       max_theta = max_fov %>% unname(),
       max_theta_px = max_fov_px)
}
