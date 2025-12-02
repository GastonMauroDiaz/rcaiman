#' Tune sky-sampling parameters
#'
#' @description
#' Evaluate combinations of two discrete parameters governing sky-point
#' extraction: the number of equal-area sky segments (`n_cells`) and the minimum
#' allowed distance to the binarized canopy (`dist_to_black`).
#' For each parameter pair, the function computes a composite performance metric
#' combining (i) angular uniformity relative to a reference segmentation and
#' (ii) local sky variability (CV) derived from the original image.
#' The pair with the lowest metric is selected.
#'
#' @param r numeric [terra::SpatRaster-class] of one layer. Typically the blue
#'   band of a canopy photograph. Digital numbers must be linearly related to
#'   radiance. See [read_caim_raw()] for details.
#' @param n_cells_seq numeric vector. Sequence of integer values to evaluate for
#'   `n_cells`.
#' @param dist_to_black_seq numeric vector. Sequence of integer values to
#'   evaluate for `dist_to_black`.
#'
#' @inheritParams compute_canopy_openness
#' @inheritParams apply_by_direction
#'
#' @details
#' Each combination (`n_cells`, `dist_to_black`) is evaluated by:
#' \enumerate{
#'   \item generating an equal-area segmentation with `n_cells`;
#'   \item extracting one representative sky point per cell, subject to the
#'     minimum distance to pixels classified as non-sky;
#'   \item computing two indicators:
#'     (a) median local sky CV,
#'     (b) angular uniformity (KL divergence) relative to a fixed reference
#'         equal-area segmentation.
#' }
#'
#' The final metric is the mean of the two indicators, expressed on comparable
#' relative scales.
#'
#'
#' @return List with three components:
#' \describe{
#'   \item{`n_cells`}{Optimal value within `n_cells_seq`.}
#'   \item{`dist_to_black`}{Optimal value within `dist_to_black_seq`.}
#'   \item{`metric`}{Minimum composite metric.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' seg <- equalarea_segmentation(z, a, 30)
#' bin <- binarize_by_region(r, seg, "thr_twocorner") & select_sky_region(z, 0, 88)
#'
#' params <- tune_sky_sampling(r, z, a, bin, seq(100, 800, 50), 1:10)
#'
#' seg <- equalarea_segmentation(z, a, params$n_cells)
#' sky_points <- extract_sky_points(r, bin, seg, params$dist_to_black)
#' display_caim(caim$Blue, sky_points = sky_points)
#' }
tune_sky_sampling <- function(r, z, a, bin, n_cells_seq, dist_to_black_seq,
                    parallel = TRUE,
                    cores = NULL,
                    logical = TRUE,
                    leave_free = 1) {

  .single <- function(r, seg, ref_seg, bin, cv, dist_to_black) {
    sky_points <- extract_sky_points(r, bin, seg,
                                     dist_to_black = dist_to_black,
                                     method = "per_cell")

    .cv <- median(extract_dn(cv, sky_points, use_window = FALSE)[,3])
    .kl <- assess_sampling_uniformity(sky_points, ref_seg)
    mean(.cv, .kl)
  }

  sky_points <- fibonacci_points(z, a, 1)
  i <- extract_dn(bin, sky_points, use_window = FALSE)[,3]
  ref_seg <- equalarea_segmentation(z, a, sum(i))

  cv <- terra::focal(r, 3, sd) / terra::focal(r, 3, mean)

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }

  if (parallel) {
    row_vals <- terra::rowFromCell(r, terra::cells(r))
    col_vals <- terra::colFromCell(r, terra::cells(r))
    r_vals <- terra::values(r)
    z_vals <- terra::values(z)
    a_vals <- terra::values(a)
    bin_vals <- terra::values(bin)
    cv_vals <- terra::values(cv)
    ref_seg_vals <- terra::values(ref_seg)
  }

  if (parallel) {
    if (length(n_cells_seq) > 1) {
      .multi_equalarea_segmentation <- function(row_vals, col_vals,
                                                z_vals, a_vals, n_cells) {
        # Rebuilt rasters
        z <- terra::rast(nrows = max(row_vals, na.rm = TRUE),
                         ncols = max(col_vals, na.rm = TRUE))
        terra::ext(z) <- terra::ext(0, ncol(z), 0, nrow(z))
        # https://spatialreference.org/ref/sr-org/7589/
        terra::crs(z) <- "epsg:7589"
        terra::values(z) <- z_vals
        names(z) <- "Zenith image"

        a <- rast(z)
        terra::values(a) <- a_vals
        names(a) <- "Azimuth image"

        equalarea_segmentation(z, a, n_cells) %>% terra::values()
      }
      seg_vals_list <- .with_cluster(cores, {
        foreach(x = n_cells_seq) %dopar% {
          .multi_equalarea_segmentation(row_vals, col_vals, z_vals, a_vals, x)
        }
      })
    } else {
      seg <- equalarea_segmentation(z, a, n_cells_seq)
      seg_vals_list <-  list(terra::values(seg))
    }
    names(seg_vals_list) <- n_cells_seq
  } else {
    seg_list <- lapply(n_cells_seq, function(x) equalarea_segmentation(z, a, x))
    names(seg_list) <- n_cells_seq
  }


  if (parallel) {
    if (length(n_cells_seq) > 1) {
      .multi_equalarea_segmentation <- function(row_vals, col_vals,
                                                z_vals, a_vals, n_cells) {
        # Rebuilt rasters
        z <- terra::rast(nrows = max(row_vals, na.rm = TRUE),
                         ncols = max(col_vals, na.rm = TRUE))
        terra::ext(z) <- terra::ext(0, ncol(z), 0, nrow(z))
        # https://spatialreference.org/ref/sr-org/7589/
        terra::crs(z) <- "epsg:7589"
        terra::values(z) <- z_vals
        names(z) <- "Zenith image"

        a <- rast(z)
        terra::values(a) <- a_vals
        names(a) <- "Azimuth image"

        equalarea_segmentation(z, a, n_cells) %>% terra::values()
      }
      seg_vals_list <- .with_cluster(cores, {
        foreach(x = n_cells_seq) %dopar% {
          .multi_equalarea_segmentation(row_vals, col_vals, z_vals, a_vals, x)
        }
      })
    } else {
      seg <- equalarea_segmentation(z, a, n_cells_seq)
      seg_vals_list <-  list(terra::values(seg))
    }
    names(seg_vals_list) <- n_cells_seq
  } else {
    seg_list <- lapply(n_cells_seq, function(x) equalarea_segmentation(z, a, x))
    names(seg_list) <- n_cells_seq
  }


  params <- expand.grid(n_cells_seq, dist_to_black_seq)
  names(params) <- c("n_cells", "dist_to_black")

  if (parallel) {
    .multi <- function(row_vals, col_vals, r_vals, seg_vals, ref_seg_vals,
                       bin_vals, cv_vals, dist_to_black) {
      # Rebuilt rasters
      r <- terra::rast(nrows = max(row_vals, na.rm = TRUE),
                       ncols = max(col_vals, na.rm = TRUE))
      terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
      # https://spatialreference.org/ref/sr-org/7589/
      terra::crs(r) <- "epsg:7589"
      terra::values(r) <- r_vals

      cv <- bin <- seg <- ref_seg <- rast(r)

      terra::values(seg) <- seg_vals
      names(seg) <- "Sky segments of equal area"
      terra::values(ref_seg) <- ref_seg_vals
      names(ref_seg) <- "Sky segments of equal area"
      terra::values(bin) <- bin_vals
      bin <- binarize_with_thr(bin, 0.5)
      terra::values(cv) <- cv_vals

      .single(r, seg, ref_seg, bin, cv, dist_to_black)
    }
    metrics <- .with_cluster(cores, {
      foreach(i = 1:nrow(params)) %dopar% {
        .multi(row_vals, col_vals, r_vals,
               seg_vals_list[[as.character(params[i, 1])]], ref_seg_vals,
               bin_vals, cv_vals, params[i, 2])
      }
    })
  } else {
    metrics <- lapply(1:nrow(params),
                       function(i) {
                         .single(r, seg_list[[as.character(params[i, 1])]],
                                 bin, cv, params[i, 2])
                         }
                       )
  }
  metrics <- unlist(metrics)

  values <- params[which.min(metrics), ]
  list(n_cells = values$n_cells,
       dist_to_black = values$dist_to_black,
       metric = metrics[which.min(metrics)])
}
