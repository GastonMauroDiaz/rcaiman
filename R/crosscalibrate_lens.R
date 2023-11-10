#' Cross-calibrate lens
#'
#' Cross-calibrate lens
#'
#' Read the help page of [calibrate_lens()] for understanding the theory
#' being this function.
#'
#' This function is intended to be used when a camera calibrated with a method
#' of higher accuracy than the one proposed in [calibrate_lens()] is
#' available or there is a main camera to which all other devices should be
#' adjusted.
#'
#' It requires two photographs taken from the exact same location with the
#' calibrated and uncalibrated camera. This means that the lens entrance pupils
#' should match and the optical axes should be aligned.
#'
#' Points should be digitized in tandem with ImageJ and saved as CSV files.
#'
#' @param path_to_csv_uncal,path_to_csv_cal Character vector of length one. Path
#'   to a CSV file created with the
#'   [point
#'   selection tool of ‘ImageJ’
#'   software](https://imagej.net/ij/docs/guide/146-19.html#sec:Multi-point-Tool)
#'   (*cal* and *uncal* stand for calibrated and uncalibrated, respectively).
#' @param zenith_colrow_uncal,zenith_colrow_cal Numeric vector of length two.
#'   Raster coordinates of the zenith. See
#'   [calc_zenith_colrow()] (*cal* and *uncal* stand
#'   for calibrated and uncalibrated, respectively).
#' @param diameter_cal Numeric vector of length one. Diameter in pixels of the
#'   image taken with the calibrated camera.
#' @param lens_coef numeric
#' @inheritParams calibrate_lens
#'
#' @return An object of class *list* with named elements. *ds* is the dataset
#'   used to fit the model, *model* is the fitted model (class `lm`, see
#'   [stats::lm()]), *horizon_radius* is the radius at 90º, *lens_coef* is a
#'   numeric vector of length equal to the `degree` argument containing the
#'   polynomial model coefficients for predicting relative radius
#'   (`coefficients(model)/horizon_radius`).
#'
#' @export
#' @family Lens Functions
crosscalibrate_lens <- function(path_to_csv_uncal,
                           path_to_csv_cal,
                           zenith_colrow_uncal,
                           zenith_colrow_cal,
                           diameter_cal,
                           lens_coef,
                           degree = 3) {

  csv_uncal <- utils::read.csv(path_to_csv_uncal)
  csv_uncal <- cbind(csv_uncal$X, csv_uncal$Y)
  csv_cal <- utils::read.csv(path_to_csv_cal)
  csv_cal <- cbind(csv_cal$X, csv_cal$Y)

  ## center in (0,0)
  csv_uncal[, 1] <- csv_uncal[, 1] - zenith_colrow_uncal[1]
  csv_uncal[, 2] <- csv_uncal[, 2] - zenith_colrow_uncal[2]
  csv_cal[, 1] <- csv_cal[, 1] - zenith_colrow_cal[1]
  csv_cal[, 2] <- csv_cal[, 2] - zenith_colrow_cal[2]

  csv_uncal <- pracma::cart2pol(csv_uncal)
  csv_cal <- pracma::cart2pol(csv_cal)
  px <- csv_uncal[, 2]

  angle <- seq(0, 90,  length.out = 90)
  R <- calc_relative_radius(angle, lens_coef)
  inv_fun <- splinefun(R * diameter_cal/2, .degree2radian(angle))

  theta <- inv_fun(csv_cal[, 2])

  fit <- lm(px ~ poly(theta, degree, raw = TRUE) - 1)
  horizon_radius <- stats::predict(fit, data.frame(theta = pi / 2)) %>% unname()

  list(ds = data.frame(theta, px),
       model = fit,
       horizon_radius = horizon_radius,
       lens_coef = (coefficients(fit) / horizon_radius) %>% unname())
}
