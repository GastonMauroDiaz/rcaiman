#' Calculate relative radius
#'
#' Calculate the relative radius given a zenith angle and lens function. This is
#' known as the projection function.
#'
#' @param angle Numeric vector. Zenith angles in degrees.
#' @inheritParams zenith_image
#'
#' @export
#'
calc_relative_radius <- function(angle, lens_coef) {

  angle <- .degree2radian(angle)

  coef_order <- cbind(lens_coef, seq(1, length(lens_coef)))
  for (i in 1:length(lens_coef)) {
    if (i == 1) {
      ma <- coef_order[i, 1] * angle^coef_order[i, 2]
    } else {
      ma <- rbind(ma, coef_order[i, 1] * angle^coef_order[i, 2])
    }
  }

  if (length(lens_coef) == 1) {
    relative_radius <- ma
  } else {
    relative_radius <- apply(ma, 2, sum)
  }
  unname(relative_radius)
}
