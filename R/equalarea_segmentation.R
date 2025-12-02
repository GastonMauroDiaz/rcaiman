#' Segment a hemisphere into equal-area cells in spherical space
#'
#' @description
#' Segment a hemispherical view into \code{n_cells} equal-area cells in
#' spherical space.
#'
#' @details
#' Segmentation begins by creating `n_cells` initial rings through uniform
#' sampling in `cos(z)` space. The core logic of the algorithm is grouping `k`
#' initial rings into a larger ring and then subdividing that ring azimuthally
#' into `k` cells. This procedure ensures equal area while allowing the final
#' segments to adopt compact shapes. In this way, the larger rings become the
#' final rings, of which there will be `n_rings = round(sqrt(n_cells))`.
#'
#' The final ring containing the zenith has `k = 1`, the ring touching the
#' horizon has `k = n_rings`, and the intermediate `k` values must be
#' determined. To impose a compactness structure across zenith angles, the
#' method uses the equisolid-angle projection:
#'
#' \deqn{
#' r = 2 \, f \, \sin\!\left(\frac{\theta}{2}\right)
#' }
#'
#' where \eqn{r} is the radial distance, \eqn{f} is the focal length, and
#' \eqn{\theta} is the zenith angle. Based on this transformation, the
#' expression `round(2 * sin(theta / 2) * n_rings * x)` is used to compute `k`,
#' with `x` being a scalar (estimated with `optim()`) that adjusts the sum of
#' `k` values to exactly match `n_cells`.
#'
#' Because the true central zenith angles of the final rings are not known in
#' advance (because they depend on `k`) the procedure is iterated.
#' It starts by assuming rings of equal angular width to obtain provisional
#' central zenith angles, computes `k`, updates the ring boundaries from the
#' resulting grouping, recomputes the central zenith angles, and repeats until
#' convergence.
#'
#' @param z [terra::SpatRaster] single-layer raster with zenith angle (degrees).
#' @param a [terra::SpatRaster] single-layer raster with azimuth angle (degrees).
#' @param n_cells numeric scalar. Target number of cells.
#'
#' @return Single-layer [terra::SpatRaster] with integer values from 1 to \code{n_cells}.
#'
#' @seealso [skygrid_segmentation()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(500, lens())
#' a <- azimuth_image(500, lens())
#' seg <- equalarea_segmentation(z, a, n_cells = 512)
#' plot(seg)
#' }
equalarea_segmentation <- function(z, a, n_cells) {
  .check_r_z_a_m(NULL, z, a)
  .check_vector(n_cells, "integerish", 1, sign = "positive")

  initial_ring_breaks <- acos(seq(0,1, length.out = n_cells + 1))
  initial_ring_breaks <- initial_ring_breaks %>%  .radian2degree() %>% rev()
  initial_ring_breaks <- initial_ring_breaks[-1]
  n_rings <- round(sqrt(n_cells))

  .transform.c_theta <- function(x) {
    c_theta <- .degree2radian(c_theta)
    ############################################################################
    # equisolid-angle transformation
    # Schneider, D., Schwalbe, E., & Maas, H.-G. (2009). Validation of geometric
    # models for fisheye lenses. ISPRS Journal of Photogrammetry and Remote
    # Sensing, 64(3), Article 3. https://doi.org/10.1016/j.isprsjprs.2009.01.001
    v <- 2 * sin(c_theta / 2)
    ############################################################################
    round(v * n_rings * x)
  }

  # calc central theta assuming equal-angle rings
  delta_theta <- 90 / n_rings
  c_theta <- seq(0 + delta_theta / 2, 90 - delta_theta / 2, by = delta_theta)

  .fn <- function(x) {
    n_cells_per_ring <- .transform.c_theta(x)
    abs(n_cells - sum(n_cells_per_ring))
  }
  delta <- 1
  repeat {
    # Define number of cells per ring
    x <- tryCatch(optim(1, .fn, method = "SANN")$par, error = function(e) 1)
    n_cells_per_ring <- .transform.c_theta(x)
    rest <- n_cells - sum(n_cells_per_ring)
    n_cells_per_ring[n_rings] <- n_cells_per_ring[n_rings] + rest

    # final rings
    ring_breaks <- initial_ring_breaks[cumsum(n_cells_per_ring)]
    ring_breaks <- c(0, ring_breaks)

    #new central theta
    pre_c_theta <- c_theta
    c_theta <- apply(data.frame(ring_breaks[-length(ring_breaks)],
                                    ring_breaks[-1]), 1, mean)
    # to continue or not to continue
    delta <- sqrt(sum(c_theta - pre_c_theta)^2)
    if (delta < 1) break
  }
  if (rest > 0) message(paste0("The last ring has ", rest, " cells more than `n_rings`."))

  # compose the raster
  rcl <- data.frame(ring_breaks[-length(ring_breaks)], ring_breaks[-1], 1:n_rings)
  rings <- terra::classify(z, rcl)

  gain <- cumsum(n_cells_per_ring)
  gain <- c(0, gain)
  seg <- lapply(seq_along(n_cells_per_ring),
                   function(i) {
                      m <- rings == i
                      if (n_cells_per_ring[i] == 1) {
                        sectors <- m
                      } else {
                        sectors <- terra::classify(a, n_cells_per_ring[i])
                        sectors <- as.numeric(sectors) + 1
                      }
                      sectors[!m] <- NA
                      sectors + gain[i]
                     }
                )
  seg <- rast(seg)
  seg <- terra::cover(seg)

  names(seg) <- "Sky segments of equal area"
  seg
}
