#' Lens database
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
#' modeled with a straight line. Following
#' \href{https://www.schleppi.ch/patrick/hemisfer/}{Hemisfer software}, this
#' package uses a polynomial curve to model lens distortion. A third-order
#' polynomial is sufficient in most cases
#' \insertCite{Frazer2001;textual}{rcaiman}.
#'
#' Eventually, this will be a large database, but only the following lenses are
#' available at the moment:
#'
#' \itemize{ \item \strong{equidistant}: standard equidistant projection
#' \insertCite{Schneider2009}{rcaiman}.
#'
#' \item \strong{Nikon_FCE9}: Nikon FC-E9 auxiliary lens
#' \insertCite{Diaz2018}{rcaiman}
#'
#' \item \strong{Nikkor_10.5_mm}: AF DX Fisheye-Nikkor 10.5mm f/2.8G ED
#' \insertCite{Pekin2009}{rcaiman}
#'
#' \item \strong{Olloclip}: Auxiliary lens. Unpublished
#'
#' }
#'
#' @param type Character vector of length one. The name of the lens.
#' @param max_fov Logical vector of length one. Use \code{TRUE} to return the
#'   maximum field of view in degrees.
#'
#' @export
#'
#' @return If \code{max_fov} is set to \code{TRUE}, it returns a numeric vector
#'   of length one, which is the lens maximum field of view in degrees.
#'   Otherwise, it returns a numeric vector with the coefficients of the lens
#'   function.
#'
#' @family Lens Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' lens("Nikon_FCE9")
#' lens("Nikon_FCE9", max_fov = TRUE)
lens <- function(type = "equidistant", max_fov = FALSE) {
  if (max_fov) index <- 2 else index <- 1

  type <- trimws(type)

  switch(type,
    equidistant = list(2 / pi, 180)[[index]],
    Nikon_FCE9 = list(c(0.6427, 0.0346, -0.024491), 190)[[index]],
    Nikkor_10.5_mm = list(c(0.71553, 0.01146, -0.03928), 165)[[index]],
    Soligor_fisheye = list(c(0.6427, 0.0346, -0.024491), 180)[[index]],
    Olloclip = list(c(0.7836131, 0.1511661, -0.1558095), 170)[[index]]
  )
}
