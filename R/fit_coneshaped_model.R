#' Fit cone-shaped model
#'
#' Generate the digital numbers of the whole sky through statistical modelling.
#'
#' An explanation of this function can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}. However, although the model is
#' the same, this implementation is more flexible thank to
#' \code{\link{extract_sky_points}} and \code{\link{extract_zenith_dn}}.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018}{rcaiman}.
#'
#' @param sky_points The data.frame returned by \code{\link{extract_zenith_dn}},
#'   or a data.frame with the same structure and names.
#' @param use_azimuth_angle Logical vector of length one. If \code{TRUE},
#'   Equation 4 from \insertCite{Diaz2018;textual}{rcaiman} is used: \eqn{sDN =
#'   a + b \cdot \theta + c  \cdot \theta^2 + d  \cdot sin(\phi) + e  \cdot
#'   cos(\phi)}, where \eqn{sDN} is sky digital number, \eqn{a,b,c,d} and
#'   \eqn{e} are coefficients, \eqn{\theta} is zenith angle, and \eqn{\phi} is
#'   azimuth angle. If \code{FALSE}, a simplified version based on
#'   \insertCite{Wagner2001;textual}{rcaiman} is used: \eqn{sDN = a + b \cdot
#'   \theta + c  \cdot \theta^2}.
#'
#' @return A list of two objects, one of class \code{function} and the other of
#'   class \code{lm} (see \code{\link[stats]{lm}}). If the fitting fails, it
#'   returns \code{NULL}. The function requires two arguments, azimuth and
#'   zenith in degrees, and returns relative luminance.
#' @export
#'
#' @family MBLT functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels(blue, z, a)
#' sky_points <- extract_sky_points(blue, bin, g)
#' zenith_dn <- extract_zenith_dn(blue, z, a, sky_points)
#' rl_cs_fun <- fit_coneshaped_model(zenith_dn$sky_points)
#' sky_cs <- rl_cs_fun$rl_cs_fun(a, z) * zenith_dn$zenith_dn
#' persp(sky_cs, theta = 90, phi = 0) #a flipped rounded cone!
#' }
fit_coneshaped_model <- function(sky_points,
                                 use_azimuth_angle = TRUE) {
  stopifnot(class(sky_points) == "data.frame")
  stopifnot(length(use_azimuth_angle) == 1)

  Blue <- sky_points$rl
  Zenith <- sky_points$a
  Azimuth <- sky_points$z

  if (length(Blue) > 30) {
    if (use_azimuth_angle) {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE) +
        sin(Azimuth * pi / 180) + cos(Azimuth * pi / 180))

      # only to avoid note from check, code is OK without this line.
      a <- b <- d <- e <- NA

      skyFun <- function(z, azimuth) {
        x <- coefficients(model)
        x[is.na(x)] <- 0
        for (i in 1:5) assign(letters[i], x[i])
        a + b * z + c * z^2 +
          d * sin(azimuth * pi / 180) + e * cos(azimuth * pi / 180)
      }
    } else {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE))
      skyFun <- function(z, azimuth) {
        x <- coefficients(model)
        x[is.na(x)] <- 0
        for (i in 1:5) assign(letters[i], x[i])
        a + b * z + c * z^2
      }
    }
    return(list(rl_cs_fun = skyFun, model = model))
  } else {
    return(NULL)
  }
}
