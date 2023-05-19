crosscalibrate_lens <- function(path_to_csv_uncal,
                           path_to_csv_cal,
                           diameter,
                           lens_coef,
                           degree = 3) {

  csv_uncal <- utils::read.csv(path_to_csv_uncal)
  csv_uncal <- cbind(csv_uncal$X, csv_uncal$Y)
  csv_cal <- utils::read.csv(path_to_csv_cal)
  csv_cal <- cbind(csv_cal$X, csv_cal$Y)

  ## center in (0,0)
  csv_uncal[, 1] <- csv_uncal[, 1] - csv_uncal[1, 1]
  csv_uncal[, 2] <- csv_uncal[, 2] - csv_uncal[1, 2]
  csv_cal[, 1] <- csv_cal[, 1] - csv_cal[1, 1]
  csv_cal[, 2] <- csv_cal[, 2] - csv_cal[1, 2]

  csv_uncal <- pracma::cart2pol(csv_uncal)
  px <- csv_uncal[, 2]
  max_fov_px <- max(px)

  angle <- seq(0, 90,  length.out = 90)
  R <- calc_relative_radius(angle, lens_coef)
  inv_fun <- splinefun(R * diameter/2, .degree2radian(angle))

  theta <- inv_fun(csv_cal[, 2])

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
