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
#' @param sky_points The *data.frame* returned by [extract_rel_radiance()] or a
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
#'   relative radiance
#' @export
#'
#' @seealso [thr_mblt()]
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' com <- compute_complementary_gradients(caim)
#' chroma <- max(com$blue_yellow, com$cyan_red)
#' bin <- apply_thr(chroma, thr_isodata(chroma[!is.na(chroma)]))
#' bin <- bin & apply_thr(com$blue_yellow, -0.2)
#'
#' g <- sky_grid_segmentation(z, a, 10, first_ring_different = TRUE)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_black = 3)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' rr <- extract_rel_radiance(r, z, a, sky_points)
#'
#' model <- fit_coneshaped_model(rr$sky_points)
#' summary(model$model)
#' sky_cs <- model$fun(z, a) * rr$zenith_dn
#' plot(sky_cs)
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

  Blue <- sky_points$rr
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
    warning("Insufficient number of points to attempt model fitting")
    return(NULL)
  }
}
