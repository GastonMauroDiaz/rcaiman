#' Lens database
#'
#' Database of lens projection functions and field of views.
#'
#' Eventually, this will be a large database, but only next lens
#' are available for the moment:
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
#' }
#'
#' @param type Character vector of length one. The name of the lens, see
#'   details.
#' @param max_fov Logical. Use \code{TRUE} to return the maximum field of view
#'   in degrees.
#'
#' @export
#'
#' @family Lens functions
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#'   \insertRef{Pekin2009}{rcaiman}
#'
#'   \insertRef{Schneider2009}{rcaiman}
#'
#' @examples lens("equidistant")
lens <- function(type = "equidistant", max_fov = FALSE) {
  if (max_fov) index <- 2 else index <- 1

  type <- trimws(type)

  switch(type,
    equidistant = list(2 / pi, 180)[[index]],
    Nikon_FCE9 = list(c(0.6427, 0.0346, -0.024491), 190)[[index]],
    Nikkor_10.5_mm = list(c(0.71553, 0.01146, -0.03928), 165)[[index]],
    Soligor_fisheye = list(c(0.6427, 0.0346, -0.024491), 180)[[index]],
    Olloclip = list(c(1.06065, -0.49054, 0.14044), 165)[[index]]
  )
}
