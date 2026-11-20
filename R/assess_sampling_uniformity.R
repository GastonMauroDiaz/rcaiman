#' Assess sampling uniformity
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Quantify the distribution of sky samples within a segmented field of view
#' using Kullback--Leibler divergence and angular-area coverage.
#'
#' @details
#' The reference domain is the field of view of the oriented acquisition,
#' represented by positive labels in `seg` with positive, available `w` values.
#' An optional `m` further restricts this domain. For sky reconstruction, use
#' a segmentation of external directions, such as [beckers_segmentation()] or
#' [healpix_segmentation()], which exclude directions below the horizon.
#'
#' The sky gaps supply samples within this domain; they do not define it.
#' Keep `seg`, `w`, and `m` fixed when comparing sampling alternatives. In
#' particular, do not use each candidate sky mask as `m`: regions obscured by
#' vegetation still belong to the domain to be reconstructed.
#'
#' For each segment, sum the pixel weights within the reference domain to obtain
#' its relative angular area \eqn{A_i}. Let \eqn{c_i} be its number of samples.
#' The empirical and reference proportions are
#' \deqn{p_i = c_i / \sum_j c_j, \qquad q_i = A_i / \sum_j A_j.}
#' Each sampling row contributes one count, including repeated indices. Weights
#' describe the reference areas; they do not weight the sample counts.
#'
#' Kullback--Leibler divergence is
#' \deqn{\mathrm{KL}(p\|q) = \sum_{i:p_i>0} p_i \log(p_i/q_i).}
#' Coverage is the proportion of reference area in segments receiving samples:
#' \deqn{\mathrm{coverage} = \sum_{i:c_i>0} q_i.}
#' Thus, coverage depends on segmentation resolution; it is not the area of
#' the sampled pixels themselves. When reference segments have equal areas,
#' these expressions reduce to equal expected counts and a fraction of cells.
#'
#' @note
#' Equal-area constructors define equal complete cells, but the FOV, horizon,
#' or a deliberate restriction can cut those cells. Summing `w` accounts for
#' their represented portions, including rasterization. Unobserved cells are
#' not part of the reference. A represented segment without samples remains
#' part of it. Compare alternatives using the same domain and segmentation;
#' the measures do not guarantee the accuracy of a subsequent sky model.
#'
#' @param seg single-layer [terra::SpatRaster-class] of nonnegative integer-like
#'   labels. Positive labels identify segments; zero and `NA` are excluded.
#'   IDs need not be consecutive. No constructor-specific attributes are needed.
#' @param w single-layer numeric [terra::SpatRaster-class] on the same grid as
#'   `seg`, with nonnegative finite solid-angle weights, for example from
#'   [area_weight_image()]. Relative weights are sufficient. Missing and zero
#'   weights are excluded. The reference must contain positive weight.
#' @param m optional single-layer logical, NA-free [terra::SpatRaster-class]
#'   on the same grid as `seg`. `TRUE` retains directions in the reference
#'   domain. With `NULL`, all positive labels with positive weights are used.
#'   This is a restriction of the domain to reconstruct, not a sky-gap mask.
#' @param sampling_indices `data.frame` with positive integer-like `col` and
#'   `row` columns on the grid of `seg`. Every row must select a pixel within
#'   the reference domain. Other columns are ignored. Use [rem_unselected()]
#'   to filter indices beforehand when needed.
#'
#' @return List with `kl`, the Kullback--Leibler divergence from angular-area
#'   proportions, and `coverage`, the fraction of reference area in segments
#'   receiving at least one sample. With no samples, `kl` is `NA` and coverage
#'   is zero.
#' @export
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
#' w <- area_weight_image(g, geometry)
#' seg <- beckers_segmentation(d, 20)
#' # Alternatively: seg <- healpix_segmentation(d, nside = 2)
#' sampling_indices <- expand.grid(col = seq(1, 80, 4), row = seq(1, 60, 4))
#' sampling_indices <- rem_unselected(sampling_indices, seg > 0)
#' assess_sampling_uniformity(sampling_indices, seg, w)
#'
#' # Restrict the reconstruction domain, keeping it fixed across comparisons.
#' m <- !is.na(d$zenith) & d$zenith <= 70
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' assess_sampling_uniformity(sampling_indices, seg, w, m)
assess_sampling_uniformity <- function(sampling_indices, seg, w, m = NULL) {
  reference <- .sampling_area_reference(seg, w, m)
  .assess_sampling_reference(sampling_indices, reference)
}

# Build once per optimization so every candidate uses the same domain.
.sampling_area_reference <- function(seg, w, m = NULL) {
  .assert_single_layer(seg)
  labels <- .segment_values(seg)
  .assert_area_weight_image(w)
  .assert_same_geom(seg, w)
  weights <- terra::values(w, mat = FALSE)
  selected <- !is.na(labels) & labels > 0 & !is.na(weights) & weights > 0
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(seg, m)
    selected <- selected & terra::values(m, mat = FALSE)
  }
  if (!any(selected)) {
    stop("The reference domain must contain positive weights and positive cell IDs.",
         call. = FALSE)
  }
  weights <- weights[selected] / max(weights[selected])
  areas <- rowsum(weights, labels[selected], reorder = TRUE)
  ids <- as.numeric(rownames(areas))
  q <- as.numeric(areas) / sum(areas)
  labels[!selected] <- NA_real_
  domain <- terra::setValues(terra::rast(seg), labels)
  names(domain) <- "segment"
  list(seg = domain, ids = ids, q = q)
}

.sampling_reference_ids <- function(sampling_indices, reference) {
  .check_columns(sampling_indices, required_columns = c("col", "row"))
  indices <- sampling_indices[, c("col", "row"), drop = FALSE]
  ids <- extract_dn(reference$seg, indices, use_window = FALSE)[[3]]
  if (anyNA(ids)) {
    stop("Sampling indices must select positive cell IDs within the reference domain.",
         call. = FALSE)
  }
  match(ids, reference$ids)
}

.assess_sampling_reference <- function(sampling_indices, reference) {
  ids <- .sampling_reference_ids(sampling_indices, reference)
  if (!length(ids)) return(list(kl = NA_real_, coverage = 0))
  p <- tabulate(ids, nbins = length(reference$q)) / length(ids)
  occupied <- p > 0
  list(kl = sum(p[occupied] * log(p[occupied] / reference$q[occupied])),
       coverage = sum(reference$q[occupied]))
}
