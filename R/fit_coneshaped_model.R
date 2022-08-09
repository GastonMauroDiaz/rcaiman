#' Fit cone-shaped model
#'
#' Statistical modeling to predict the digital numbers from spherical
#' coordinates.
#'
#' An explanation of this model can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018;textual}{rcaiman} in addition to this package.
#'
#' @param sky_points The \emph{data.frame} returned by \code{\link{extract_rl}}
#'   or a \emph{data.frame} with same structure and names.
#' @param use_azimuth_angle Logical vector of length one. If \code{TRUE}, the
#'   Equation 4 from \insertCite{Diaz2018;textual}{rcaiman}) is used: \eqn{sDN =
#'   a + b \cdot \theta + c  \cdot \theta^2 + d  \cdot sin(\phi) + e  \cdot
#'   cos(\phi)}, where \eqn{sDN} is sky digital number, \eqn{a,b,c,d} and
#'   \eqn{e} are coefficients, \eqn{\theta} is zenith angle, and \eqn{\phi} is
#'   azimuth angle. If \code{FALSE}, the next simplified version based on
#'   \insertCite{Wagner2001;textual}{rcaiman} is used: \eqn{sDN = a + b \cdot
#'   \theta + c  \cdot \theta^2}.
#'
#' @return A list of two objects, one of class \code{function} and the other of
#'   class \code{lm} (see \code{\link[stats]{lm}}). If the fitting fails, it
#'   returns \code{NULL}. The function requires two arguments--zenith and
#'   azimuth in degrees--to return relative luminance.
#' @export
#'
#' @family Sky Reconstruction Functions
#' @seealso \code{\link{thr_image}}
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels(r, z, a)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_rl(r, z, a, sky_points, NULL)
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' sky_cs <- model$fun(z, a)
#' persp(sky_cs, theta = 90, phi = 0) #a flipped rounded cone!
#' }
fit_coneshaped_model <- function(sky_points,
                                 use_azimuth_angle = TRUE) {
  stopifnot(class(sky_points) == "data.frame")
  stopifnot(length(use_azimuth_angle) == 1)

  Blue <- sky_points$rl
  Zenith <- sky_points$z
  Azimuth <- sky_points$a

  if (length(Blue) > 30) {
    if (use_azimuth_angle) {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE) +
        sin(Azimuth * pi / 180) + cos(Azimuth * pi / 180))

      # only to avoid note from check, code is OK without this line.
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
