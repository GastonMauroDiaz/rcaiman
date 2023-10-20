.gaussian2d <- function(x, y, target_a, target_b, sigma) {
  stats::dnorm(x, target_a, sigma) * stats::dnorm(y, target_b, sigma)
}

.get_gaussian_2d_parameters <- function(target_color, sigma) {
  if (!is(target_color, "LAB")) {
    if (!is(target_color, "sRGB")) {
      target_color <- as(target_color, "sRGB")
    } else {
      target_color <- as(target_color, "LAB")
    }
  }
  ma <- colorspace::coords(target_color)
  target_a <- ma[, 2]
  target_b <- ma[, 3]
  if (is.null(sigma)) sigma <- sqrt(target_a^2 + target_b^2) %>% unname()
  c(target_a, target_b, sigma)
}

#' Compute the membership to a target color
#'
#' This function was first presented in \insertCite{Diaz2015;textual}{rcaiman}.
#' It computes the degree of membership to a color using two Gaussian membership
#' functions and the axes _a*_ and _b*_ from the
#' _CIE L*a*b*_ color space. The lightness dimension is not
#' considered in the calculations.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2015;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`).
#'
#' @inheritParams expand_noncircular
#' @param target_color [color-class].
#' @param sigma Numeric vector of length one. Use `NULL` (default) to estimate
#'   it automatically as the euclidean distance between `target_color` and grey
#'   in the _CIE L*a*b*_ color space.
#'
#'
#' @return It returns an object from the class [SpatRaster-class]. First layer
#'   is the membership to the target color. Second layer is the membership to
#'   grey. Both memberships are calculated with same `sigma`.
#'
#' @export
#' @family Pre-processing Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' plotRGB((caim/(2^16-2))*255)
#' sky_blue_sample <- crop_caim(caim, c(327, 239), 41, 89)
#' plotRGB((sky_blue_sample/2^16)*255)
#' sky_blue <- apply(sky_blue_sample[], 2, median) %>% normalize(., 0, 2^16) %>%
#'   as.numeric() %>%
#'   matrix(., ncol = 3) %>%
#'   sRGB()
#' caim <- normalize(caim)
#' mem <- membership_to_color(caim, sky_blue)
#' plot(mem)
#' }
membership_to_color <- function(caim, target_color, sigma = NULL) {
  .is_class_from_colorspace(target_color)
  stopifnot(class(caim) == "SpatRaster")
  stopifnot(all(names(caim) == c("Red", "Green", "Blue")))
  .was_normalized(caim, "caim")
  stopifnot(names(caim) == c("Red", "Green", "Blue"))
  if (!is.null(sigma)) stopifnot(length(sigma) == 1)

  color <- colorspace::sRGB(terra::values(caim))
  color <- as(color, "LAB")

  p <- .get_gaussian_2d_parameters(target_color, sigma)
  max_z <- .gaussian2d(p[1], p[2], p[1], p[2], p[3])
  x <- colorspace::coords(color)
  mem_to_color <- .gaussian2d(x[, 2], x[, 3], p[1], p[2], p[3]) / max_z

  target_color <- colorspace::sRGB(matrix(c(0, 0, 0), ncol = 3))
  sigma <- p[3]
  p <- .get_gaussian_2d_parameters(target_color, sigma)
  max_z <- .gaussian2d(p[1], p[2], p[1], p[2], sigma)
  mem_to_grey <- .gaussian2d(x[, 2], x[, 3], p[1], p[2], sigma) / max_z

  r <- terra::subset(caim, 1:2)
  r$Red <- mem_to_color
  r$Green <- mem_to_grey
  names(r) <- c("membership_to_target_color", "membership_to_grey")
  r
}
