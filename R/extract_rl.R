#' Extract relative luminance
#'
#' Extract the luminance relative to the zenith digital number
#'
#' The search for near-zenith points starts in the region  ranged between `0`
#' and `z_thr`. If the number of near-zenith points is less than `no_of_points`,
#' the region increases by steps of `2` degrees of zenith angle till the
#' required number of points is reached.
#'
#' @inheritParams ootb_mblt
#' @param sky_points An object of class *data.frame*. The output of
#'   [extract_sky_points()]. As an alternative, both
#'   [ImageJ](https://imagej.net/ij/) and HSP software package
#'   \insertCite{Lang2013}{rcaiman} can be used to manually digitize points. See
#'   [extract_dn()] and [read_manual_input()] for details.
#' @param no_of_points Numeric vector of length one. The number of near-zenith
#'   points required for the estimation of the zenith DN.
#' @param z_thr Numeric vector on length one. The starting maximum zenith angle
#'   used to search for near-zenith points.
#' @param use_window Logical vector of length one. If `TRUE`, a \eqn{3 \times 3}
#'   window will be used to extract the digital number from `r`.
#' @param min_spherical_dist Numeric vector of length one. This parameter
#'   filters out points based on their proximity on a spherical surface.
#'   Essentially, this is the spherical counterpart of the `min_raster_dist`
#'   argument from [extract_sky_points()]. The distance is expressed in degrees.
#'
#' @return A list of three objects, *zenith_dn* and *max_zenith_angle* from the
#'   class *numeric*, and *sky_points* of the class
#'   *data.frame*; *zenith_dn* is the estimated zenith digital number,
#'   *max_zenith_angle* is the maximum zenith angle reached in the search
#'   for near-zenith sky points, and *sky_points* is the input argument
#'   `sky_points` with the additional columns: *a*, *z*,
#'   *dn*, and *rl*, which stand for azimuth and zenith angle in
#'   degrees, digital number, and relative luminance, respectively. If `NULL` is
#'   provided as `no_of_points`, then *zenith_dn* is forced to one and,
#'   therefore, *dn* and *rl* will be identical.
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
#' rl <- extract_rl(caim$Blue, z, a, sky_points, 1, min_spherical_dist = 10)
#' points(rl$sky_points$col, nrow(caim) - rl$sky_points$row, col = 3, pch = 0)
#' }
extract_rl <- function(r, z, a, sky_points,
                       no_of_points = 3,
                       z_thr = 5,
                       use_window = TRUE,
                       min_spherical_dist = NULL) {

  if (!is.null(min_spherical_dist)) {
    stopifnot(length(min_spherical_dist) == 1)
    stopifnot(is.numeric(min_spherical_dist))
  }

  .dist_spherical <- function(zenith_azimuth) {
    n <- nrow(zenith_azimuth)
    m <- matrix(0, n, n)

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        d <- .calc_spherical_distance(zenith_azimuth[i, 1],
                                      zenith_azimuth[i, 2],
                                      zenith_azimuth[j, 1],
                                      zenith_azimuth[j, 2], radians = FALSE)
        m[i, j] <- d
        m[j, i] <- d
      }
    }
    dimnames(m) <- list(rownames(zenith_azimuth), rownames(zenith_azimuth))
    m
  }
  #below function was taken from extract_sky_points() and modified
  .filter <- function(ds, thr) {
    d <- .dist_spherical(ds[, c("z", "a")])
    indices <- c()
    i <- 0
    while (i < nrow(d)) {
      i <- i + 1
      indices <- c(indices, row.names(d)[i]) #include the point itself (p)
      x <- names(d[i, d[i,] <= thr])
      if (!is.null(x)) {
        # this exclude from future search all the points near p,
        # including itself
        rows2crop <- (1:nrow(d))[match(x, rownames(d))]
        cols2crop <- (1:ncol(d))[match(x, colnames(d))]
        d <- d[-rows2crop, -cols2crop]
      }
      if (is.vector(d)) d <- matrix(d)
    }
    indices
  }

  stopifnot(is.data.frame(sky_points))
  stopifnot(ncol(sky_points) == 2)
  .check_if_r_z_and_a_are_ok(r, z, a)

  cells <- terra::cellFromRowCol(a, sky_points$row, sky_points$col)
  sky_points$a <- a[cells][,]
  sky_points$z <- z[cells][,]

  if(any(is.na(sky_points$z))) {
      stop(paste0("This problem arises from having sky points touching ",
                  "the horizon. It generally solves masking out the region ",
                  "near the horizon using 'bin <- bin & mask_hs(z, 0, 80)' ",
                  "as a preprocessing step.")
           )
  }

  if (!is.null(min_spherical_dist)) {
    i <- sor_filter(sky_points, r,
                    k = 3,
                    rmax = min_spherical_dist,
                    thr = 0,
                    cutoff_side = "left")
    sky_points <- sky_points[i, ]
    sky_points <- sky_points[.filter(sky_points, min_spherical_dist %>%
                                       .degree2radian()), ]
    cells <- terra::cellFromRowCol(a, sky_points$row, sky_points$col)
  }

  if (use_window) {
    xy <-  terra::xyFromCell(r, cells)
    r_smooth <- terra::focal(r, 3, "mean")
    sky_points$dn <-  terra::extract(r_smooth, xy)[,]
    if(any(is.na(sky_points$dn))) {
      warning("The kernel created NA values.")
    }
  } else {
    sky_points$dn <-  r[cells][,]
  }

  if (is.null(no_of_points)) {
    zenith_dn <- 1
  } else {
    .is_whole(no_of_points)
    stopifnot(length(no_of_points) == 1)
    zenith_dn <- c()
    unlock <- TRUE
    while (unlock) {
      zenith_dn <- sky_points$dn[sky_points$z < z_thr]
      z_thr <- z_thr + 2
      unlock <- length(zenith_dn) < no_of_points
      if (z_thr >= 90) unlock <- FALSE
    }
    if ((z_thr - 2) > 20) warning(paste0("Zenith DN were estimated from ",
                                         "sky pixels from zenith angle up to ",
                                         z_thr))
    zenith_dn <- mean(zenith_dn)
  }
  sky_points$rl <- sky_points$dn / zenith_dn
  list(zenith_dn = zenith_dn, max_zenith_angle = z_thr, sky_points = sky_points)
}
