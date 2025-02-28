#' Calibrate lens
#'
#' Calibrate a fisheye lens
#'
#' Fisheye lenses have a wide field of view and the same distortion in all
#' directions running orthogonally to the optical axis. The latter property
#' allows fitting a precise mathematical relationship between distances to the
#' zenith on the image space and zenith angles on the hemispherical space
#' (assuming upward-looking hemispherical photography with the optical axis
#' vertically aligned).
#'
#' The method
#' outlined here, known as the simple method, is explained in details in
#' \insertCite{Diaz2024;textual}{rcaiman}. Next explanation might serve mostly
#' as a handbook.
#'
#' ## Step-by-step guide for producing a CSV file to feed this function
#'
#' ### Materials
#'
#' * this package and [ImageJ](https://imagej.net/ij/download.html)
#' * camera and lens
#' * tripod
#' * standard yoga mat
#' * table at least as wide as the yoga mat width
#' * twenty two push pins of different colors
#' * one print of this [sheet](https://osf.io/tudzc) (A1 size,
#' almost like a research poster).
#' * scissors
#' * some patience
#'
#' ### Instructions
#'
#' Cut the sheet by the dashed line. Place the yoga mat extended on top of the
#' table. Place the sheet on top of the yoga mat. Align the dashed line with the
#' yoga mat border closest to you. Place push pins on each cross. If you are
#' gentle, the yoga mat will allow you to do that without damaging the table. Of
#' course, other materials could be used to obtain the same result, such as
#' cardboard, foam, nails, etc.
#'
#' ![](calibrationBoard.jpg "Calibration board")
#'
#' Place the camera on the tripod. Align its optical axis with the table while
#' looking for getting an image showing the overlapping of the three pairs of
#' push pins, as instructed in the print. In order to take care of the line of
#' pins at 90º relative to the optical axis, it may be of help to use the naked
#' eye to align the entrance pupil of the lens with the pins. The alignment of
#' the push pins only guarantees the position of the lens entrance pupil, the
#' leveling should be cheeked with an instrument, and the alignment between the
#' optical axis and the radius of the zenith push pin should be taken into
#' account. In practice, the latter is achieved by aligning the camera body with
#' the orthogonal frame made by the quarter circle.
#'
#' Take a photo and transfer it to the computer, open it with ImageJ, and use
#' the [point selection
#' tool](https://imagej.net/ij/docs/guide/146-19.html#sec:Multi-point-Tool)
#' to digitize the push pins, starting from the zenith push pin and not skipping
#' any shown push pin. End with an additional point where the image meets the
#' surrounding black (or the last pixel in case there is not blackness because
#' it is not a circular hemispherical image. There is no need to follow the line
#' formed by the push pins). Then, use the dropdown menu Analyze>Measure to open
#' the window Results. To obtain the CSV, use File>Save As...
#'
#' ![](pushpinsImageJ.jpg "Points digitization with ImageJ")
#'
#' Use [test_lens_coef()] to test if coefficients are OK.
#'
#'
#' @note
#'
#' If we imagine the fisheye image as an analog clock, it is possible to
#' calibrate 3 o'clock by attaching the camera to the tripod in landscape mode
#' while leaving the quarter-circle at the lens's right side. To calibrate 9
#' o'clock, it will be necessary to rotate the camera to put the quarter-circle
#' at the lens's left side. To calibrate 12 and 6 o'clock, it will be necessary
#' to do the same but with the camera in portrait mode. If several directions
#' are sampled with this procedure, a character vector of length greater than
#' one in which each element is a path to a CSV files could be provided as the
#' `path_to_csv` argument.
#'
#' @param path_to_csv Character vector. Path to a CSV file created with the
#'   [point selection tool of ‘ImageJ’
#'   software](https://imagej.net/ij/docs/guide/146-19.html#sec:Multi-point-Tool).
#' @param degree Numeric vector of length one. Polynomial model degree.
#'
#' @return An object of class *list* with named elements. *ds* is the dataset
#'   used to fit the model, *model* is the fitted model (class `lm`, see
#'   [stats::lm()]), *horizon_radius* is the radius at 90º, *lens_coef* is a
#'   numeric vector of length equal to the `degree` argument containing the
#'   polynomial model coefficients for predicting relative radius
#'   (`coefficients(model)/horizon_radius`),
#'   *zenith_colrow* are the raster coordinates of the zenith push pin,
#'   *max_theta* is the maximum zenith angle in degrees, and *max_theta_px* is
#'   the distance in pixels between the zenith and the maximum zenith angle in
#'   pixels units.
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
#' coefficients(calibration$model)
#' calibration$lens_coef %>% signif(3)
#' calibration$horizon_radius
#'
#' \dontrun{
#' test_lens_coef(calibration$lens_coef) #MacOS and Windows tend to differ here
#' test_lens_coef(c(0.628, 0.0399, -0.0217))
#' }
#'
#' .fp <- function(theta, lens_coef) {
#'   x <- lens_coef[1:5]
#'   x[is.na(x)] <- 0
#'   for (i in 1:5) assign(letters[i], x[i])
#'   a * theta + b * theta^2 + c * theta^3 + d * theta^4 + e * theta^5
#' }
#'
#' plot(calibration$ds)
#' theta <- seq(0, pi/2, pi/180)
#' lines(theta, .fp(theta, coefficients(calibration$model)))
calibrate_lens <- function(path_to_csv, degree = 3) {
  .fun <- function(path_to_csv) {
    csv <- utils::read.csv(path_to_csv)
    csv <- cbind(csv$X, csv$Y)
    zenith_colrow <- csv[1,]

    # center in (0,0)
    csv[, 1] <- csv[, 1] - csv[1, 1]
    csv[, 2] <- csv[, 2] - csv[1, 2]

    csv <- pracma::cart2pol(csv)
    px <- csv[, 2]
    max_theta_px <- px[length(px)]
    # remove the data point that is only useful for max_fox_px calculation
    px <- px[-length(px)]

    theta <- .degree2radian(seq(0, 90, 5))
    theta <- theta[1:length(px)] # because some lens would "see" less pushpins
    ds <- data.frame(theta, px)
    list(ds = ds,
         zenith_colrow = zenith_colrow,
         max_theta_px = max_theta_px)
  }
  if (length(path_to_csv) > 1) {
    l <- Map(.fun, path_to_csv)
    ds <- l[[1]]$ds
    for (i in 2:length(path_to_csv)) {
      ds <- rbind(ds, l[[i]]$ds)
    }
    zenith_colrow <- Map(function(x) x$zenith_colrow, l) %>% unlist() %>%
      matrix(., ncol = 2, byrow = TRUE) %>% apply(., 2, stats::median)
    max_theta_px <- Map(function(x) x$max_theta_px, l) %>% unlist() %>%
      stats::median()
    l <- list(ds = ds,
              zenith_colrow = zenith_colrow,
              max_theta_px = max_theta_px)
  } else {
    l <- .fun(path_to_csv)
  }
  theta <- l$ds$theta
  px <- l$ds$px

  fit <- lm(theta ~ poly(px, degree, raw = TRUE) - 1)
  lens_coef_caneye <- unname(coefficients(fit))
  max_theta <- stats::predict(fit, data.frame(px = l$max_theta_px)) %>%
    .radian2degree()

  fit <- lm(px ~ poly(theta, degree, raw = TRUE) - 1)
  horizon_radius <- stats::predict(fit, data.frame(theta = pi / 2)) %>% unname()

  list(ds = l$ds,
       model = fit,
       horizon_radius = horizon_radius,
       lens_coef = (coefficients(fit) / horizon_radius) %>% unname(),
       zenith_colrow = l$zenith_colrow,
       max_theta = max_theta %>% unname(),
       max_theta_px = l$max_theta_px)
}
