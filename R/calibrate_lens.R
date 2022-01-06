#' Calibrate lens
#'
#' Calibrate a fisheye lens. This type of lens has wide field of view and a
#' consistent azimuthal distortion, so a precise mathematical relation can be
#' fit between the distance to the zenith on the image space and the zenith
#' angle on the hemispherical space
#'
#' If you cannot find the coefficient of your lens on the literature, you may
#' want to try the solution offered here. It requires, in addition to this
#' package and the open-source
#' \href{https://imagej.nih.gov/ij/download.html}{ImageJ software package}, the
#' following materials:
#'
#' \itemize{ \item camera and lens
#'
#' \item tripod
#'
#' \item standard yoga mat
#'
#' \item table larger than the yoga mat
#'
#' \item twenty two push pins of different colors
#'
#' \item scissors
#'
#' \item One print of this \href{https://osf.io/tudzc/download}{sheet} (A1 size,
#' almost like a poster).
#'
#' }
#'
#' Cut the sheet by the dashed line. Place the yoga mat extended on top of the
#' table. Place the sheet on top of the yoga mat. Align the dashed line with the
#' yoga mat border closest to you, and place push pins on each cross. If you are
#' gentle, the yoga mat will allows you to do that without damaging the table.
#' Of course, other materials could be used to obtain the same result, such as
#' cardboard, foam, nails, etc.
#'
#' Place the camera on the tripod, align its optical axis with the table while
#' looking for getting an image showing the overlapping of the three pairs of
#' push pins as instructed in the print. Take a photograph and check if it looks
#' more or less like \href{https://osf.io/tudzc/download}{this one}.
#'
#' Transfer the photograph to the computer, open it with ImageJ, and use the
#' \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#sec:Multi-point-Tool}{point
#' selection tool} to digitize the push pins, starting from the zenith push pin
#' and not skipping any showed push pin. This method was inspired by the
#' calibration board from \insertCite{Clark1988;textual}{rcaiman}.
#'
#' As a tip, use \code{\link{test_lens_coef}} to test if coefficient are OK. If
#' not, try moving the last points a little bit. Put the last one a few pixels
#' farther from the zenith is usually enough.
#'
#' @param path_to_csv Character vector of length one. Path to a CSV file created
#'   with the
#'   \href{https://imagej.nih.gov/ij/docs/guide/146-19.html#sec:Multi-point-Tool}{point
#'    selection tool of ‘ImageJ’ software}.
#' @param degree Numeric vector of length one. Polynomial model degree.
#'
#' @return An object of class list with named elements. ‘lens_coef’ stands for
#'   lens coefficients, ‘max_theta’ for maximum zenith angle in degrees, and
#'   ‘max_theta_px’ for distance in pixels between the zenith and the maximum
#'   zenith angle in pixels units.
#'
#' @family Lens functions
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
#' calibration$max_thera_px
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
       max_theta = max_fov,
       max_theta_px = max_fov_px)
}
