#' Extract relative luminance
#'
#' Extract the luminance relative to the zenith digital number.
#'
#' To estimate the zenith digital number (zenith DN), the search for near-zenith
#' points starts in the region  region ranged between \code{0} and \code{z_thr}.
#' If the number of near zenith points is less than \code{no_of_points}, the
#' region increase by \code{2} degrees of zenith angle till the required number
#' of points is reached
#'
#' @inheritParams ootb_mblt
#' @param sky_points An object of class \emph{data.frame}. The result of a call
#'   to \code{\link{extract_sky_points}}.
#' @param no_of_points Numeric vector on length one. The number of near-zenith
#'   points required for the estimation of the zenith DN.
#' @param z_thr Numeric vector on length one. The starting maximum zenith angle
#'   used to search for near zenith points.
#' @param use_window Logical vector of length one. If \code{TRUE}, a \eqn{3
#'   \times 3} window will be used to extract the sky digital number from
#'   \code{r}.
#'
#' @return A list of three objects, \emph{zenith_dn} and \emph{max_zenith_angle}
#'   from the class numeric, and \emph{sky_points} from the class
#'   \emph{data.frame}; \emph{zenith_dn} is the estimated zenith digital number,
#'   \emph{max_zenith_angle} is the maximum zenith angle reached in the search
#'   for near-zenith sky points, and \emph{sky_points} is the input argument
#'   \code{sky_points} with the additional columns: \emph{a}, \emph{z},
#'   \emph{dn}, and \emph{rl}, which stand for azimuth and zenith angle in
#'   degrees, digital number, and relative luminance, respectively. If
#'   \code{NULL} is provided as \code{no_of_points}, then \emph{zenith_dn} is
#'   forced to one and \emph{dn}, and \emph{rl} are equals.
#' @export
#'
#' @family Tools Functions
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
#' rl <- extract_rl(r, z, a, sky_points, 1)
#' }
extract_rl <- function(r, z, a, sky_points,
                       no_of_points = 20,
                       z_thr = 2,
                       use_window = TRUE) {
  stopifnot(is.data.frame(sky_points))
  .check_if_r_z_and_a_are_ok(r, z, a)

  cells <- terra::cellFromRowCol(a, sky_points$row, sky_points$col)
  sky_points$a <- a[cells][,]
  sky_points$z <- z[cells][,]

  if(any(is.na(sky_points$z))) {
      stop("Check for any white pixel outside the circle.")
    }

  if (use_window) {
    xy <-  terra::xyFromCell(r, cells)
    r_smooth <- terra::focal(r, 3, "mean")
    sky_points$dn <-  terra::extract(r_smooth, xy)[,]
    if(any(is.na(sky_points$dn))) {
      warning("The kernel created NA values near the horizon.")
    }
  } else {
    sky_points$dn <-  r[cells][,]
  }

  if (is.null(no_of_points)) {
    zenith_dn <- 1
  } else {
    .is_integerish(no_of_points)
    stopifnot(length(no_of_points) == 1)
    zenith_dn <- c()
    unlock <- TRUE
    while (unlock) {
      zenith_dn <- sky_points$dn[sky_points$z < z_thr]
      z_thr <- z_thr + 2
      unlock <- length(zenith_dn) < no_of_points
      if (z_thr >= 90) unlock <- FALSE
    }
    if ((z_thr - 2) > 20) warning(paste0("Zenith DN were estimated from pure ",
                                         "sky pixels from zenith angle up to ",
                                         z_thr))
    zenith_dn <- mean(zenith_dn)
  }

  sky_points$rl <- sky_points$dn / zenith_dn
  list(zenith_dn = zenith_dn, max_zenith_angle = z_thr, sky_points = sky_points)
}
