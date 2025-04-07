#' Extract relative radiance
#'
#' Extract the radiance relative to the zenith radiance
#'
#' Interpolate the `no_of_points` closer to the zenith using IDW with p = 2.
#'
#' @inheritParams extract_dn
#' @inheritParams ootb_mblt
#' @param no_of_points Numeric vector of length one. The number of near-zenith
#'   points required for the estimation of the zenith DN.
#' @param use_window Logical vector of length one. If `TRUE`, a \eqn{3 \times 3}
#'   window will be used to extract the digital number from `r`.
#'
#' @note As an alternative, both [ImageJ](https://imagej.net/ij/) and HSP
#' software package \insertCite{Lang2013}{rcaiman} can be used to manually
#' digitize points. See [extract_dn()] and [read_manual_input()] for details.
#'
#' @return A list of three objects, *zenith_dn* and *max_zenith_angle* from the
#'   class *numeric*, and *sky_points* of the class
#'   *data.frame*; *zenith_dn* is the estimated zenith digital number,
#'   *max_zenith_angle* is the maximum zenith angle reached in the search
#'   for near-zenith sky points, and *sky_points* is the input argument
#'   `sky_points` with the additional columns: *a*, *z*,
#'   *dn*, and *rr*, which stand for azimuth and zenith angle in
#'   degrees, digital number, and relative radiance, respectively. If `NULL` is
#'   provided as `no_of_points`, then *zenith_dn* is forced to one and,
#'   therefore, *dn* and *rr* will be identical.
#' @export
#'
#' @family Tool Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # See fit_cie_sky_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' rr <- extract_rel_radiance(caim$Blue, z, a, sky_points, 1)
#' points(rr$sky_points$col, nrow(caim) - rr$sky_points$row, col = 3, pch = 0)
#' }
extract_rel_radiance <- function(r, z, a, sky_points,
                       no_of_points = 3,
                       use_window = TRUE) {

  stopifnot(is.data.frame(sky_points))
  stopifnot(ncol(sky_points) == 2)
  if (is.null(no_of_points)) {
    stopifnot(nrow(sky_points) >= 1)
  } else {
    stopifnot(nrow(sky_points) >= no_of_points)
  }
  .check_if_r_z_and_a_are_ok(r, z, a)

  # Extract spherical coordinates
  cells <- terra::cellFromRowCol(a, sky_points$row, sky_points$col)
  sky_points$a <- a[cells][,]
  sky_points$z <- z[cells][,]
  if(any(is.na(sky_points$z))) {
      stop(paste0("This problem arises from having sky points touching ",
                  "the horizon. It generally solves masking out the region ",
                  "near the horizon using ",
                  "'bin <- bin & select_sky_vault_region(z, 0, 80)' ",
                  "as a preprocessing step.")
           )
  }

  # Extract values from the image with the points
  if (use_window) {
    xy <-  terra::xyFromCell(r, cells)
    r_smoothed <- terra::focal(r, 3, "mean")
    sky_points$dn <-  terra::extract(r_smoothed, xy)[,]
    if(any(is.na(sky_points$dn))) {
      warning("The kernel created NA values.")
    }
  } else {
    sky_points$dn <-  r[cells][,]
  }

  # Estimate zenith DN
  if (is.null(no_of_points)) {
    zenith_dn <- 1
  } else {
    spherical_distance <- .calc_spherical_distance(
                                              sky_points$z %>% .degree2radian(),
                                              sky_points$a %>% .degree2radian(),
                                              0,
                                              0,
                                              radians = TRUE)
    k <- no_of_points
    sorted_indices <- order(spherical_distance)
    w <- spherical_distance[sorted_indices][2:(k + 1)]
    w <- 1 / w^2
    u <- sky_points[sorted_indices[2:(k + 1)], "dn"]
    zenith_dn <- sum(u * (w / sum(w)))
  }

  # Calculate relative radiance
  sky_points$rr <- sky_points$dn / zenith_dn


  list(zenith_dn = zenith_dn,
       sky_points = sky_points)
}
