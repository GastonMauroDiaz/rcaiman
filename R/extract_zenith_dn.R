#' Extract zenith digital number
#'
#' The search for near zenith points starts in the region delimited by
#' \code{z_thr}. If the number of near zenith points is less than
#' \code{no_of_points}, the region increase by 2 degrees of zenith angle til the
#' required number of points is reached
#'
#' @inheritParams ootb_mblt
#' @param sky_points An object of class data.frame. The result of a call to
#'   \code{\link{extract_sky_points}}.
#' @param no_of_points An object of class numeric. The number of near zenith
#'   points required for the estimation of the zenith dn.
#' @param z_thr An object of class numeric. The starting maximum zenith angle
#'   used to search for near zenith points.
#' @param use_window Logical vector of length one. If \code{TRUE}, a 3 by 3
#'   window will be used to extract the sky digital number from \code{r}.
#'
#' @return A list of three objects, ‘zenith_dn’ and ‘max_zenith_angle’ from the
#'   class numeric, and ‘sky_points’ from the class data.frame; ‘zenith_dn’ is
#'   the estimated zenith digital number, ‘max_zenith_angle’ is the maximum
#'   azimuth reached in the search for near zenith sky points, and ‘sky_points’
#'   is the input argument ‘sky_points’ with the additional columns: \emph{a},
#'   \emph{z}, and \emph{dn}, and \emph{rl}, which stand for azimuth and zenith
#'   angle in degrees, digital number, and relative luminance, respectively.
#' @export
#'
#' @family MBLT functions
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
#' zenith_dn <- extract_zenith_dn(blue, z, a, sky_points, 1)
#' }
extract_zenith_dn <- function(r, z, a, sky_points,
                              no_of_points = 20,
                              z_thr = 2,
                              use_window = TRUE) {
  stopifnot(is.data.frame(sky_points))
  .check_if_r_z_and_a_are_ok(r, z, a)

  cells <- terra::cellFromRowCol(a, sky_points$row, sky_points$col)
  sky_points$a <- a[cells][,]
  sky_points$z <- z[cells][,]

  if (use_window) {
    xy <-  terra::xyFromCell(r, cells)
    r_smooth <- terra::focal(r, 3, "mean")
    sky_points$dn <-  terra::extract(r_smooth, xy)[,]
  } else {
    sky_points$dn <-  r[cells][,]
  }

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
  sky_points$rl <- sky_points$dn / zenith_dn
  list(zenith_dn = zenith_dn, max_zenith_angle = z_thr, sky_points = sky_points)
}
