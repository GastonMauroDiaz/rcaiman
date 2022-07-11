#' Lens database
#'
#' Database of lens projection functions and field of views.
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
#' }
#'
#' @param type Character vector of length one. The name of the lens, see
#'   details.
#' @param max_fov Logical. Use \code{TRUE} to return the maximum field of view
#'   in degrees.
#'
#' @export
#'
#' @return If \code{max_fov} is set to \code{TRUE}, it returns a numeric vector
#'   of length one, which is the lens maximum field of view in degrees.
#'   Otherwise, it returns a numeric vector with the coefficient of the lens
#'   function.
#'
#' @family Lens functions
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
