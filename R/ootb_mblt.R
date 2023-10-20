#' Out-of-the-box model-based local thresholding
#'
#' Out-of-the-box version of the model-based local thresholding (MBLT) algorithm
#'
#' This function is a hard-coded version of a MBLT pipeline. The MBLT approach
#' proposes a linear relationship between background value and optimal threshold
#' value. This function uses statistical models for sky reconstruction that are
#' able to explain smooth changes in sky brightness, so this function works best
#' under clear skies or overcast conditions. After the reconstruction, local
#' thresholds are linearly predicted from sky brightness values (see
#' [thr_mblt()]).
#'
#' As a high-level summary, the function starts producing a working binarized
#' image and ends with a refined binarized image.
#'
#' The pipeline combines these main functions [extract_sky_points_simple()] or
#' [extract_sky_points()], [fit_coneshaped_model()], and [fit_trend_surface()].
#' The code can be easily inspected by calling `ootb_mblt` without parenthesis.
#' Advanced users can use that code as a template.
#'
#' The MBLT algorithm was first presented in
#' \insertCite{Diaz2018;textual}{rcaiman}. The version presented here differs
#' from the original in the following main aspects:
#'
#' * The original version used a global thresholding method to provide sky
#' points to the cone-shaped model. This one uses [extract_sky_points_simple()].
#' Nevertheless, a binarized image can be provided through the `bin` argument
#' ([find_sky_pixels()], [regional_thresholding()], [ootb_obia(), and
#' [ootb_mblt()]] are recommended ways to obtain a working binarized image). If
#' a binarized image is provided, [extract_sky_points()] is used.
#' * `intercept` and `slope` are automatically obtained using data from sky
#' points and a linear model for accuracy evaluation after
#' \insertCite{Pineiro2008;textual}{rcaiman}. This approach handles inaccuracies
#' in background reconstruction (see [thr_mblt()] for additional details).
#' * This version does not use asynchronous acquisition under the open sky, as
#' the original method did. The cone-shaped model ([fit_coneshaped_model()]) run
#' without a filling source and the cone-shaped sky is used as filling source
#' for trend surface fitting ([fit_trend_surface()]).
#'
#' This function searches for black objects against a light background. When
#' regular canopy hemispherical images are provided as input, the algorithm will
#' find dark canopy elements against a bright sky almost everywhere in the
#' picture and, therefore, the result will fit user expectations. However, if a
#' hemispherical photograph taken under the open sky is provided, this algorithm
#' would be still searching black objects against a light background, so the
#' darker portions of the sky will be taken as objects, i.e., canopy. As a
#' consequence, this will not fit users expectations since they are looking for
#' the classes *Gap* and *No-gap*, no matter if one of those are not in the
#' picture itself. This kind of error could happen with photographs of open
#' forests for the same working principle.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018;textual}{rcaiman} in addition to this package
#' (`citation("rcaiman"`).
#'
#' @param r [SpatRaster-class]. A normalized greyscale image. Typically, the
#'   blue channel extracted from a canopy photograph. Please see [read_caim()]
#'   and [normalize()].
#' @param z [SpatRaster-class] built with [zenith_image()].
#' @param a [SpatRaster-class] built with [azimuth_image()].
#' @param bin [SpatRaster-class]. This should be a preliminary binarization of
#'   `r` useful for masking pixels that are very likely pure sky pixels.
#' @param fix_cs_sky Logical vector of length one. If it is `TRUE`,
#'   [fix_reconstructed_sky()] is used to fix the cone-shaped sky.
#' @param w Numeric vector of length one. Weighting parameter from
#'   \insertCite{Diaz2018;textual}{rcaiman}'s Equation 1. See [thr_mblt()]
#'
#' @note
#'
#' If `NULL` is provided as the `w` argument, the weight is calculated as the
#' coefficient of determination (\eqn{R^2}) of linear model for accuracy
#' evaluation \insertCite{Pineiro2008}{rcaiman}.
#'
#' @export
#' @family Binarization Functions
#'
#' @return Object from class list containing the binarized image (named
#'   *bin*) and the reconstructed skies (named *sky_cs* and
#'   *sky_s*).
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' r <- correct_vignetting(r, z, c(0.0638, -0.101)) %>% normalize()
#' bin <- find_sky_pixels(r, z, a)
#' bin <- ootb_mblt(r, z, a, bin)
#' plot(bin$bin)
#'
#' }
ootb_mblt <- function(r, z, a, bin = NULL, fix_cs_sky = FALSE, w = 0.5) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  .fast <- function() {
    sky_points <- extract_sky_points_simple(r, z, a)
    sky_points <- extract_rl(r, z, a, sky_points, NULL, use_window = FALSE)
    fit_coneshaped_model(sky_points$sky_points)
  }

  g <- sky_grid_segmentation(z, a, 10)
  if (!is.null(bin)) {
    try(sky_points <- extract_sky_points(r, bin, g), silent = TRUE)
    if (exists("sky_points")) {
      sky_points <- extract_rl(r, z, a, sky_points, NULL)
      model <- fit_coneshaped_model(sky_points$sky_points)
    } else {
      model <- .fast()
    }
  } else {
    model <- .fast()
  }

  if (is.null(model)) {
    sky_cs <- z
    terra::values(sky_cs) <- median(sky_points$sky_points$dn)
    mblt <- c(0,1)
  } else {
    sky_cs <- model$fun(z, a)
    x <- predict(model$model)
    y <- predict(model$model) + model$model$residuals
    mblt <- coefficients(lm(x~y))
  }

  if (fix_cs_sky) sky_cs <- fix_reconstructed_sky(sky_cs, z, r, bin)
  sky_cs <- normalize(sky_cs, 0, 1, TRUE)
  bin <- apply_thr(r, thr_mblt(sky_cs, mblt[1]*255, mblt[2] * 0.5))
  sky_s <- fit_trend_surface(r, z, a, bin,
                             filling_source = sky_cs,
                             np = 6)$image

  sky_points <- extract_sky_points(r, bin, g)
  x <- extract_dn(sky_s, sky_points)[,3]
  y <- extract_dn(r, sky_points)[,3]
  fit <- lm(x~y)
  mblt <- coefficients(fit)
  if (is.null(w)) w <- summary(fit)$r.squared

  sky_s <- normalize(sky_s, 0, 1, TRUE)
  bin <- apply_thr(r, thr_mblt(sky_s, mblt[1]*255, mblt[2]*w))
  list(bin = bin, sky_cs = sky_cs, sky_s = sky_s)
}
