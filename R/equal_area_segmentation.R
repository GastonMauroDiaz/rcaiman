beckers_equal_area_rings <- function(n_cells, p = pi, min_k0 = 10) {

  k <- n_cells
  k_seq <- c()

  while (k >= min_k0) {
    k_seq <- c(k, k_seq)
    next_k <- round((sqrt(k) - sqrt(p))^2)

    if (next_k < min_k0 || next_k == 0) {
      break
    }
    k <- next_k
  }

  if (k_seq[1] > 1) {
    k_seq <- c(1, k_seq)
  }

  n_cells_per_ring <- diff(c(0, k_seq))

  ring_breaks_rad <- acos(1 - k_seq / n_cells)
  ring_breaks_deg <- ring_breaks_rad * 180 / pi

  list(
    n_rings = length(n_cells_per_ring),
    n_cells_per_ring = n_cells_per_ring,
    cumulative_cells = k_seq,
    ring_breaks_deg = ring_breaks_deg,
    ring_breaks_rad = ring_breaks_rad
  )
}


#' Segment a hemisphere into equal-area cells in spherical space
#'
#' @description
#' Partition a hemispherical view into equal-area cells in spherical space.
#'
#' @details
#' This function implements the algebraic recurrence method proposed by
#' \insertCite{Beckers2012;textual}{rcaiman} for partitioning a hemisphere into
#' equal-area cells. The method derives concentric rings from cumulative cell
#' counts and subdivides each ring into azimuthal sectors so that all cells have
#' the same solid angle.
#'
#' @param n_cells numeric vector of length one. Must be positive and interpreted
#'   as integer. Target number of cells.
#' @param shape_control numeric vector of length one. The aspect ratio parameter
#'   \eqn{p} from \insertCite{Beckers2012;textual}{rcaiman}. Default is
#'   \code{pi}, optimizing for square-like cells.
#' @param min_k0 numeric vector of length one. Minimum cumulative cell count
#'   before forcing the unitary zenith cell.
#'
#' @inheritParams add_directions
#'
#' @return Single-layer [terra::SpatRaster] with integer values from 1 to
#'   \code{n_cells}.
#'
#' @seealso [healpix_segmentation()]
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # effect of the projection on the apparent shape
#' z <- zenith_image(500, equal_area_coef)
#' a <- azimuth_image(z)
#' seg <- equal_area_segmentation(z, a, n_cells = 20)
#' plot(seg)
#'
#' z <- zenith_image(500, equidistant_coef)
#' a <- azimuth_image(z)
#' seg <- equal_area_segmentation(z, a, n_cells = 20)
#' plot(seg)
#'
#' z <- zenith_image(500, orthographic_coef)
#' a <- azimuth_image(z)
#' seg <- equal_area_segmentation(z, a, n_cells = 20)
#' plot(seg)
#'
#' # effect of shape_control on the actual shape
#' z <- zenith_image(500, equal_area_coef)
#' a <- azimuth_image(z)
#' seg <- equal_area_segmentation(z, a, n_cells = 20, shape_control = pi)
#' plot(seg)
#' seg <- equal_area_segmentation(z, a, n_cells = 20, shape_control = 0.5)
#' plot(seg)
#' }
equal_area_segmentation <- function(z, a, n_cells, shape_control = pi, min_k0 = 3) {
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(n_cells, "integerish", 1, sign = "positive")

  params <- beckers_equal_area_rings(n_cells, p = shape_control, min_k0 = min_k0)
  ring_breaks <- c(0, params$ring_breaks_deg)
  n_cells_per_ring <- params$n_cells_per_ring
  n_rings <- length(params$n_cells_per_ring)

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

  seg[is.na(seg)] <- 0
  names(seg) <- "Sky segments of equal area"
  attr(seg, "ring_breaks") <- ring_breaks
  seg
}
