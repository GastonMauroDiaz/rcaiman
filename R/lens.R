#' Access the lens database
#'
#' Database of lens projection functions and field of views.
#'
#' In upward-looking leveled hemispherical photography, the zenith is the center
#' of a circle whose perimeter is the horizon. This is true only if the lens
#' field of view is 180º. The relative radius is the radius of concentric
#' circles expressed as a fraction of the radius that belongs to the circle that
#' has the horizon as perimeter. The equidistant model, also called polar, is
#' the most widely used as a standard reference. Real lenses can approximate the
#' projection models, but they always have some kind of distortion. In the
#' equidistant model, the relation between zenith angle and relative radius is
#' modeled with a straight line. Following [Hemisfer
#' software](https://www.schleppi.ch/patrick/hemisfer/), this package uses a
#' polynomial curve to model lens distortion. A third-order polynomial is
#' sufficient in most cases \insertCite{Frazer2001}{rcaiman}. Equations should
#' be fitted with angles in radians.
#'
#' Eventually, this will be a large database, but only the following lenses are
#' available at the moment:
#'
#' * equidistant: standard equidistant projection
#' \insertCite{Schneider2009}{rcaiman}.
#' * Nikkor_10.5mm: AF DX Fisheye Nikkor 10.5mm f/2.8G ED
#' \insertCite{Pekin2009}{rcaiman}
#' * Nikon_FCE9: Nikon FC-E9 converter
#' \insertCite{Diaz2024}{rcaiman}
#' * Olloclip: Auxiliary lens for mobile devices made by Olloclip
#' \insertCite{Diaz2024}{rcaiman}
#' * Nikkor_8mm: AF–S Fisheye Nikkor 8–15mm f/3.5–4.5E ED
#' \insertCite{Diaz2024}{rcaiman}
#'
#' @param type Character vector of length one. The name of the lens.
#' @param max_fov Logical vector of length one. Use `TRUE` to return the maximum
#'   field of view in degrees.
#'
#' @export
#'
#' @return If `max_fov` is set to `TRUE`, it returns a numeric vector of length
#'   one, which is the lens maximum field of view in degrees. Otherwise, it
#'   returns a numeric vector with the coefficients of the lens function.
#'
#' @family Lens Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' lens("Nikon_FCE9")
#' lens("Nikon_FCE9", max_fov = TRUE)
#'
#' .fp <- function(theta, lens_coef) {
#'   x <- lens_coef[1:5]
#'   x[is.na(x)] <- 0
#'   for (i in 1:5) assign(letters[i], x[i])
#'   a * theta + b * theta^2 + c * theta^3 + d * theta^4 + e * theta^5
#' }
#'
#' theta <- seq(0, pi/2, pi/180)
#' plot(theta, .fp(theta, lens()), type = "l", lty = 2,
#'       ylab = "relative radius")
#' lines(theta, .fp(theta, lens("Nikon_FCE9")))
#'
lens <- function(type = "equidistant", max_fov = FALSE) {
  if (max_fov) index <- 2 else index <- 1

  type <- trimws(type)

  switch(type,
    equidistant = list(2 / pi, 180)[[index]],
    Nikkor_10.5mm = list(c(0.716, 0.0115, -0.0393), 165)[[index]],
    Nikon_FCE9 = list(c(0.643, 0.0346, -0.0245), 190)[[index]],
    Olloclip = list(c(0.801, 0.178, -0.179), 172)[[index]],
    Nikkor_8mm = list(c(0.689, 0.0131, -0.0295), 186)[[index]]
  )
}
