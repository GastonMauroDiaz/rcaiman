beckers_equal_area_rings <- function(n_cells, p = pi, min_k0 = 10) {

  k <- n_cells
  k_seq <- c()

  repeat {
    k_seq <- c(k, k_seq)
    if (k < min_k0) break
    next_k <- round((sqrt(k) - sqrt(p))^2)

    if (next_k < min_k0 || next_k == 0) {
      break
    }
    if (next_k >= k) {
      stop("`shape_control` and `min_k0` must give decreasing ring counts.",
           call. = FALSE)
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


#' Segment the sky using the Beckers equal-area partition
#'
#' @description
#' Partition the sky hemisphere into equal-area cells and label the observed
#' external directions on the image grid.
#'
#' @details
#' This function implements the algebraic recurrence method proposed by
#' \insertCite{Beckers2012;textual}{rcaiman} for partitioning a hemisphere into
#' equal-area cells. The method derives concentric rings from cumulative cell
#' counts and subdivides each ring into azimuthal sectors so that all cells have
#' the same solid angle. Ring boundaries and sectors are defined on the complete
#' sky hemisphere, independently of image coverage or acquisition orientation.
#'
#' Ring and sector intervals include their upper boundary; the first also
#' includes zero. Missing directions and zenith angles greater than 90 deg
#' receive zero. The unitary zenith cell does not require azimuth.
#'
#' @note
#' Equal area refers to the complete cells on the sky hemisphere, each with
#' solid angle \eqn{2\pi / n_{cells}} steradians. If `d` covers only part of the
#' sky, some cells will be only partly represented and others may be absent. The
#' observed portions therefore need not have equal areas.
#'
#' Cell IDs refer to the complete partition and, consequently, the number of
#' positive IDs present can be smaller than `n_cells`.
#'
#' The raster grid also approximates cell boundaries. Even fully represented
#' cells can have slightly different estimated solid angles. For a radial projection,
#' use [area_weight_image()] with [compute_features()] to estimate the relative
#' observed areas. [assess_sampling_uniformity()] uses these areas to account
#' for partial cells when evaluating sky sampling, for both this method and
#' [healpix_segmentation()].
#'
#' @param n_cells numeric vector of length one. Must be positive and interpreted
#'   as integer. Number of cells in the complete sky hemisphere partition.
#' @param shape_control numeric vector of length one. The aspect ratio parameter
#'   \eqn{p} from \insertCite{Beckers2012;textual}{rcaiman}. Default is
#'   \code{pi}, optimizing for square-like cells.
#' @param min_k0 numeric vector of length one. Minimum cumulative cell count
#'   before forcing the unitary zenith cell. Must be a positive integer.
#'
#' @inheritParams sky_grid_segmentation
#'
#' @return Single-layer [terra::SpatRaster-class] on the same grid as `d`,
#'   with IDs from 1 to `n_cells` and zero outside the segmentation.
#'   The `ring_breaks` attribute stores zenith boundaries in degrees.
#'   The `n_cells` attribute stores the number of cells in the complete
#'   hemisphere partition, including cells absent from the image.
#'
#' @seealso [healpix_segmentation()]
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' seg <- beckers_segmentation(d, n_cells = 20)
#' plot(seg)
beckers_segmentation <- function(d, n_cells, shape_control = pi, min_k0 = 3) {
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  .check_vector(n_cells, "integerish", 1, sign = "positive")
  .check_vector(shape_control, "numeric", 1, sign = "positive")
  .check_vector(min_k0, "integerish", 1, sign = "positive")

  params <- beckers_equal_area_rings(n_cells, p = shape_control, min_k0 = min_k0)
  ring_breaks <- c(0, params$ring_breaks_deg)
  counts <- params$n_cells_per_ring
  rings <- terra::classify(z,
    cbind(head(ring_breaks, -1), tail(ring_breaks, -1), seq_along(counts)),
    right = TRUE, include.lowest = TRUE, others = 0)
  ring_ids <- terra::values(rings, mat = FALSE)
  azimuth <- terra::values(a, mat = FALSE)
  labels <- numeric(length(ring_ids))
  offset <- 0
  for (i in seq_along(counts)) {
    selected <- !is.na(ring_ids) & ring_ids == i
    if (counts[i] == 1) {
      labels[selected] <- offset + 1
    } else {
      selected <- selected & !is.na(azimuth)
      bounds <- seq(0, 360, length.out = counts[i] + 1)
      labels[selected] <- offset + as.integer(cut(azimuth[selected], bounds,
        right = TRUE, include.lowest = TRUE))
    }
    offset <- offset + counts[i]
  }
  seg <- terra::rast(z)
  terra::values(seg) <- labels
  names(seg) <- "Beckers segmentation"
  attr(seg, "ring_breaks") <- ring_breaks
  attr(seg, "n_cells") <- n_cells
  seg
}
