#' Expand sky points
#'
#' Expand sky points using a _k_ nearest neighbors approach
#'
#' @inheritParams apply_thr
#' @inheritParams sky_grid_segmentation
#' @inheritParams extract_rel_radiance
#' @inheritParams sor_filter
#' @inheritParams interpolate_sky_points
#'
#' @returns An object of the class _data.frame_. It is the input argument
#'   `sky_points` with the following additional data:
#' \itemize{
#'   \item Rows resulting of interpolating the `sky_points` argument.
#'   \item Columns _a_, _z_, _dn_, and _initial_.
#'      \itemize{
#'        \item _a_, the azimuthal angle.
#'        \item _z_:, the zenithal angle.
#'        \item _dn_, the digital number.
#'        \item _initial_, if 'TRUE' the point is from the `sky_points` argument.
#'      }
#' }
#'
#' @family Tool Functions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' bin <- apply_thr(r, thr_isodata(r[]))
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
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points2 <- expand_sky_points(r, z, a, sky_points, k = 3, angle_width = 5)
#' plot(bin)
#' points(sky_points2$col, nrow(caim) - sky_points2$row, col = 2, pch = 10)
#' }
expand_sky_points <- function(r, z, a, sky_points, angle_width = 3,
                              k = 20, p = 2, rmax = 20) {
  .this_requires_EBImage()

  stopifnot(length(k) == 1)
  stopifnot(is.numeric(k))
  stopifnot(length(p) == 1)
  stopifnot(is.numeric(p))
  stopifnot(length(rmax) == 1)
  stopifnot(is.numeric(rmax))
  rmax <- .degree2radian(rmax)

  g <- sky_grid_segmentation(z, a, angle_width, first_ring_different = TRUE)
  g <- terra::focal(g, 3, function(x) length(unique(x)))
  bin <- g == 4
  bin <- bin & select_sky_vault_region(z, 0, 90 - angle_width * 0.9)

  bwlabels <- EBImage::bwlabel(as.array(bin))
  bwlabels <- terra::setValues(bin, bwlabels)

  sky_points2 <- extract_sky_points(r, bin, bwlabels, dist_to_black = NULL)

  sky_points2 <- extract_rel_radiance(r, z, a, sky_points2, NULL)$sky_points
  sky_points <- extract_rel_radiance(r, z, a, sky_points, NULL)$sky_points
  ds <- extract_dn(r, sky_points[, c("row", "col")])

  .calculate_dn <- function(i) {
    spherical_distance <- calc_spherical_distance(sky_points$z %>%
                                                               .degree2radian(),
                                                  sky_points$a %>%
                                                               .degree2radian(),
                                                  sky_points2[i, "z"] %>%
                                                               .degree2radian(),
                                                  sky_points2[i, "a"] %>%
                                                               .degree2radian())
    sorted_indices <- order(spherical_distance)
    w <- spherical_distance[sorted_indices][2:(k + 1)]
    m <- w <= rmax
    if (all(w <= rmax)) {
      w <- 1 / w^p
      u <- ds[sorted_indices[2:(k + 1)], 3]
      return(sum(u * (w / sum(w))))
    } else {
      return(NA)
    }
  }

  new_value <- Map(.calculate_dn, seq_len(nrow(sky_points2))) %>% unlist()

  sky_points2$dn <- new_value
  sky_points2 <- sky_points2[!is.na(new_value),]
  initial <- c(rep(TRUE, nrow(sky_points)), rep(FALSE, nrow(sky_points2)))
  sky_points <- rbind(sky_points, sky_points2)
  sky_points <- sky_points[, -ncol(sky_points)] #remove the rl column
  cbind(sky_points, initial)
}
