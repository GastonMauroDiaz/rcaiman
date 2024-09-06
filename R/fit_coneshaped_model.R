#' Fit cone-shaped model
#'
#' Statistical modeling for predicting digital numbers from spherical
#' coordinates.
#'
#' This method was presented in \insertCite{Diaz2018;textual}{rcaiman}, under the
#' heading *Estimation of the sky DN as a previous step for our method*. If you
#' use this function in your research, please cite that paper in addition to
#' this package (`citation("rcaiman"`).
#'
#' @param sky_points The *data.frame* returned by [extract_rl()] or a
#'   *data.frame* with same structure and names.
#' @param use_azimuth_angle Logical vector of length one. If `TRUE`, the
#'   Equation 4 from \insertCite{Diaz2018;textual}{rcaiman}) is used: \eqn{sDN =
#'   a + b \cdot \theta + c  \cdot \theta^2 + d  \cdot sin(\phi) + e  \cdot
#'   cos(\phi)}, where \eqn{sDN} is sky digital number, \eqn{a,b,c,d} and
#'   \eqn{e} are coefficients, \eqn{\theta} is zenith angle, and \eqn{\phi} is
#'   azimuth angle. If `FALSE`, the next simplified version based on
#'   \insertCite{Wagner2001;textual}{rcaiman} is used: \eqn{sDN = a + b \cdot
#'   \theta + c  \cdot \theta^2}.
#'
#' @return A list of two objects, one of class `function` and the other of class
#'   `lm` (see [stats::lm()]). If the fitting fails, it returns `NULL`. The
#'   function requires two arguments--zenith and azimuth in degrees--to return
#'   relative luminance.
#' @export
#'
#' @family Sky Reconstruction Functions
#' @seealso [thr_mblt()]
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
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 30), "thr_isodata")
#' bin <- bin & mask_hs(z, 0, 80)
#' sky_points <- extract_sky_points(r, bin, sky_grid_segmentation(z, a, 3))
#' sky_points <- extract_rl(r, z, a, sky_points, no_of_points = NULL)
#'
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' summary(model$model)
#' sky_cs <- model$fun(z, a)
#' plot(r/sky_cs)
#' plot(sky_cs)
#' plot(r/sky_cs > 0.5)
#'
#' z <- zenith_image(50, lens())
#' a <- azimuth_image(z)
#' sky_cs <- model$fun(z, a)
#' persp(sky_cs, theta = 90, phi = 20)
#' }
fit_coneshaped_model <- function(sky_points,
                                 use_azimuth_angle = TRUE) {
  stopifnot(class(sky_points) == "data.frame")
  stopifnot(length(use_azimuth_angle) == 1)

  Blue <- sky_points$rl
  Zenith <- sky_points$z
  Azimuth <- sky_points$a

  if (length(Blue) >= 20) {
    if (use_azimuth_angle) {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE) +
        sin(Azimuth * pi / 180) + cos(Azimuth * pi / 180))

      # Only to avoid note from check, code is OK without this line.
      a <- b <- d <- e <- NA

      skyFun <- function(zenith, azimuth) {
        x <- coefficients(model)
        x[is.na(x)] <- 0
        for (i in 1:5) assign(letters[i], x[i])
        a + b * zenith + c * zenith^2 +
          d * sin(azimuth * pi / 180) + e * cos(azimuth * pi / 180)
      }
    } else {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE))
      skyFun <- function(zenith, azimuth) {
        x <- coefficients(model)
        x[is.na(x)] <- 0
        for (i in 1:5) assign(letters[i], x[i])
        a + b * zenith + c * zenith^2
      }
    }
    return(list(fun = skyFun, model = model))
  } else {
    return(NULL)
  }
}
