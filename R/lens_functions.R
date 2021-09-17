#' Lens database
#'
#' Database of lens projection functions and field of views.
#'
#' To do: include the character needed to retrieve the lens
#' type and the reference
#'
#' @param type Character vector of length one. The name of the lens, see details.
#' @param max_fov Logical. Use TRUE to return the maximum FOV in degrees.
#'
#' @export
#'
#' @examples lens("equiangular")
lens <- function(type = "equiangular", max_fov = FALSE) {
  if (max_fov) index <- 2 else index <- 1

  type <- trimws(type)

  switch(type,
    equiangular = list(2 / pi, 180)[[index]],
    Nikon_FCE9 = list(c(0.6427, 0.0346, -0.024491), 190)[[index]],
    Nikkor_10.5_mm = list(c(0.71553, 0.01146, -0.03928), 165)[[index]],
    Soligor_fisheye = list(c(0.6427, 0.0346, -0.024491), 180)[[index]],
    Olloclip = list(c(1.06065, -0.49054, 0.14044), 165)[[index]]
  )
}



#' Test lens projection functions
#'
#' Test that lens projection function works between the 0-to-1 range.
#'
#' @inheritParams zenith_image
#'
#' @export
#'
#' @examples
#' test_lens_coef(lens("Nikon_FCE9"))
#' test_lens_coef(pi / 2)
#' test_lens_coef(c(1.06065, -0.49054, 0.14044))
test_lens_coef <- function(lens_coef) {
  testthat::test_that(
    "Test that lens projection function works between the 0-to-1 range",
    {
      testthat::expect_equal(calc_relative_radius(0, lens_coef), 0)
      testthat::expect_equal(calc_relative_radius(90, lens_coef), 1)
    }
  )
}



#' Calibrate lens
#'
#' Calibrate a fisheye lens.
#'
#' @param path_to_csv Character vector of length one. Path to a CSV file created
#'   with ImageJ.
#' @param degree Numeric vector of length one. Polynomial model degree.
#'
#' @return named list.
#' @export
#'
#' @examples
#' path <- system.file("external/Results_calibration.csv", package = "rcaiman")
#' calibration <- calibrate_lens(path)
#' calibration$lens_coef
#' calibration$max_fov
#' calibration$max_fov_px
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

  theta <- degree2radian(seq(0, 90, 5))
  theta <- theta[1:length(px)] # because some lens would "see" less pushpins

  fit <- lm(theta ~ poly(px, degree, raw = TRUE) - 1)
  lens_coef_caneye <- unname(coefficients(fit))
  max_fov <- stats::predict(fit, data.frame(px = max_fov_px)) %>%
    radian2degree()

  fit <- lm(px ~ poly(theta, degree, raw = TRUE) - 1)
  px90 <- stats::predict(fit, data.frame(theta = pi / 2))

  R <- px / px90
  fit <- lm(R ~ poly(theta, degree, raw = TRUE) - 1)

  list(lens_coef = unname(coefficients(fit)),
       lens_coef_caneye = lens_coef_caneye,
       max_fov = max_fov,
       max_fov_px = max_fov_px)
}



#' Calculate zenith raster coordinates
#'
#' @inheritParams calibrate_lens
#'
#' @export
#'
#' @examples
#' path <- system.file("external/points_over_perimeter.csv",
#'                     package = "rcaiman")
#' calc_zenith_raster_coordinates(path)
calc_zenith_raster_coordinates <- function(path_to_csv) {

  if (!requireNamespace("conicfit", quietly = TRUE)) {
    stop(paste("Package \"conicfit\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  requireNamespace("conicfit", quietly = TRUE)

  x <- utils::read.csv(path_to_csv)[, -(1:5)]

  # each tracked hole have two columns
  stopifnot(ncol(x) / 2 == round(ncol(x) / 2))

  n_holes <- ncol(x) / 2

  index <- seq(1, ncol(x), 2)

  circle <- list()
  for (i in 1:n_holes) {
    circle[[i]] <- x[, c(index[i], index[i] + 1)] %>%
      .[!is.na(.[, 1]), ] %>%
      as.matrix() %>%
      conicfit::CircleFitByKasa() %>%
      .[-3]
  }

  zenith_coordinates <- matrix(circle, ncol = 2, byrow = TRUE)
  colnames(zenith_coordinates) <- c("col", "row")
  zenith_coordinates
}


#' Calculate diameter
#'
#' Calculate the diameter of a 180º fisheye image.
#'
#' This function is useful to handle devices with field of view different than
#' 180 degrees. Given a lens projection function and data points consisting of
#' radii (pixels) and their correspondent zenith angle, it returns the radius of
#' the horizon (i.e., the radius for the zenith angle equal to 90 degrees).
#'
#' It is particularly useful when working with non-circular hemispherical
#' photography. It will help to find the diameter that a circular image would
#' have if the equipment would depict the whole hemisphere.
#'
#' @inheritParams zenith_image
#' @param radius_px Numeric vector. Distance in pixels from the zenith.
#' @param angle Numeric vector. Zenith angle in degrees.
#' @export
#'
#' @examples # Nikon D50 and Fisheye Nikkor 10.5 mm lens
#' calc_diameter(lens("Nikkor_10.5_mm"), 1202, 53)
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

  if (diameter / 2 != round(diameter / 2)) {
    diameter <- diameter + 1
  }

  diameter
}

