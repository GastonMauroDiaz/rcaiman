#' Tune sky-sampling parameters
#'
#' @description
#' Tune the parameters of [sample_sky_points()].
#'
#' @param r numeric [terra::SpatRaster-class] of one layer. Typically the blue
#'   band of a canopy photograph. Digital numbers should be linearly related to
#'   radiance.
#' @param bin_list list of logical [terra::SpatRaster-class] objects of one
#'   layer. Each element is a binarization of `r`.
#' @param n_cells_seq numeric vector. Sequence of integers to evaluate for
#'   `n_cells`.
#' @param dist_to_black_seq numeric vector. Sequence of integers to evaluate for
#'   `dist_to_black`.
#' @param w numeric vector of length one. Weight controlling the balance between
#'   coverage and accuracy (see *Details*).
#' @param write_log_in Character path (without extension) used to generate
#'   two output files summarizing the tuning run. See *Details*.
#'
#' @inheritParams compute_canopy_openness
#' @inheritParams assess_sampling_uniformity
#' @inheritParams apply_by_direction
#'
#' @details
#' This function evaluates combinations of three discrete parameters for sampling
#' sky digital numbers (\eqn{\delta_\text{sky}}) using
#' [sample_sky_points()]:
#' \itemize{
#'   \item the binarized image selected from `bin_list`, provided as argument
#'     `bin`;
#'   \item the parameter `n_cells` used to generate a segmentation with
#'     [equalarea_segmentation()], provided as argument `seg`;
#'   \item the minimum allowed distance to black pixels in the binarized image,
#'     provided as argument `dist_to_black`.
#' }
#'
#' Two criteria are computed for each parameter combination:
#' \itemize{
#'   \item *accuracy*, estimated from the local variability of
#'     \eqn{\delta_\text{sky}};
#'   \item *coverage*, the fraction of equal-area cells containing at least one
#'     sampled sky point.
#' }
#'
#' The computation of these criteria is described below, together with their
#' normalization and the construction of the final comparison metric.
#'
#' \describe{
#'  \item{*Accuracy estimation*}{
#'    True sky regions tend to show smooth radiance gradients, whereas bright
#'    canopy elements or mixed pixels introduce stronger local fluctuations.
#'    Accuracy is therefore quantified with a robust local coefficient of
#'    variation computed by [extract_cv()]. For each parameter combination, the
#'    final accuracy value is the median of these local CV estimates across all
#'    sampled sky points.
#'  }
#'  \item{*Coverage estimation*}{
#'    Coverage is computed with [assess_sampling_uniformity()], i.e., the
#'    fraction of cells in `equalarea_seg` that contain at least one sky point.
#'  }
#'  \item{*Scaling of components*}{
#'    Because the magnitudes of accuracy and coverage may differ substantially
#'    across the parameter grid, both components are rescaled to \eqn{[0,1]} using
#'    [normalize_minmax()]. This ensures that the relative contribution of each
#'    component to the optimization reflects the balance specified by \eqn{w}
#'    rather than artifacts of scale.
#'  }
#'  \item{*Coverage-Accuracy Metric (CAM)*}{
#'    After normalization, the Coverage-Accuracy Metric is defined as:
#'
#'    \deqn{
#'    \mathrm{CAM} =
#'    (1 - \mathrm{coverage}_n) \cdot w
#'    \;+\;
#'    \mathrm{accuracy}_n \cdot (1 - w) \; ,
#'    }
#'
#'    were \eqn{n} is the nth parameter combination being tested.
#'
#'    Lower CAM values indicate parameter combinations that achieve both high
#'    hemispherical coverage and low local radiance variability at the sampled sky
#'    points. The function returns the combination that minimizes CAM.
#'  }
#' }
#'
#' If `write_log_in` is not `NULL`, the function
#'   writes:
#'   * `<write_log_in>.txt`: a human-readable log containing timestamps,
#'     elapsed time (minutes), system information, the tested parameter
#'     grids, the selected combination, and its performance metrics.
#'   * `<write_log_in>_full_metrics.csv`: a tabular file listing all tested
#'     parameter combinations (`params`) together with their associated
#'     metrics.
#'   If `NULL`, no files are written.
#'
#'
#' @return List with the following numeric vectors of length one:
#' \describe{
#'   \item{`n_cells`}{Optimal value within `n_cells_seq`.}
#'   \item{`dist_to_black`}{Optimal value within `dist_to_black_seq`.}
#'   \item{`bin_list_index`}{Index of the optimal binarization in `bin_list`.}
#'   \item{`accuracy`}{Median of local CV values (MAD / median) at the
#'     sampled sky points. Lower is better.}
#'   \item{`coverage`}{Fraction of reference cells with
#'     at least one sampled sky point. Higher is better.}
#'   \item{`cam`}{Minimum Coverage-Accuracy Metric (CAM) for the
#'     selected combination.}
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
#' bin_list <- lapply(c("thr_twocorner", "thr_isodata"),
#'                    function(x) {
#'                      binarize_by_region(r, seg, x) & select_sky_region(z, 0, 88)
#'                    })
#'
#' equalarea_seg <- equalarea_segmentation(z, a, 10)
#' params <- tune_sky_sampling(r, z, a,
#'                             equalarea_seg,
#'                             bin_list,
#'                             n_cells_seq = seq(100, 500, 100),
#'                             dist_to_black_seq = 1:10,
#'                             w = 0.5,
#'                             parallel = TRUE)
#'
#' sky_points <- sample_sky_points(r,
#'                                 bin_list[[params$bin_list_index]],
#'                                 equalarea_segmentation(z, a, params$n_cells),
#'                                 params$dist_to_black)
#' display_caim(caim$Blue, sky_points = sky_points)
#' }
tune_sky_sampling <- function(r, z, a,
                              equalarea_seg,
                              bin_list,
                              n_cells_seq,
                              dist_to_black_seq,
                              w = 0.5,
                              write_log_in = NULL,
                              parallel = TRUE,
                              cores = NULL,
                              logical = TRUE,
                              leave_free = 1) {

  t0 <- Sys.time()

  .check_r_z_a_m(r, z, a, m = NULL)
  .assert_single_layer(equalarea_seg)
  .check_vector(n_cells_seq, "integerish", sign = "positive")
  .check_vector(dist_to_black_seq, "integerish", sign = "positive")
  .check_vector(w, "numeric", 1, sign = "nonnegative")
  if (w > 1) (stop("`w` cannot be greater than one."))
  .check_vector(write_log_in, "character", allow_null = TRUE)
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }

  .single <- function(r, bin, seg, equalarea_seg, dist_to_black) {
    sky_points <- sample_sky_points(r, bin, seg,
                                     dist_to_black = dist_to_black,
                                     method = "per_cell")
    acc <- stats::median(extract_cv(r, sky_points))
    cov <- assess_sampling_uniformity(sky_points, equalarea_seg)$coverage
    c(acc, cov)
  }

  params <- expand.grid(seq_along(bin_list), n_cells_seq, dist_to_black_seq)
  names(params) <- c("bin_list_index", "n_cells", "dist_to_black")

# prepare for parallel processing -----------------------------------------

  if (parallel) {
    row_vals <- terra::rowFromCell(r, terra::cells(r))
    col_vals <- terra::colFromCell(r, terra::cells(r))
    r_vals <- terra::values(r)
    z_vals <- terra::values(z)
    a_vals <- terra::values(a)
    bin_vals_list <- lapply(bin_list, function(x) terra::values(x))
    equalarea_seg_vals <- terra::values(equalarea_seg)
  }

# generate a segmentation per n_cell --------------------------------------

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
        # Only to avoid note from check, code is OK without this line.
        x <- NA

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


# calculate the metric ----------------------------------------------------

  if (parallel) {
    .multi <- function(row_vals, col_vals, r_vals, bin_vals,
                       seg_vals, equalarea_seg_vals,
                       dist_to_black) {
      # Rebuilt rasters
      r <- terra::rast(nrows = max(row_vals, na.rm = TRUE),
                       ncols = max(col_vals, na.rm = TRUE))
      terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
      # https://spatialreference.org/ref/sr-org/7589/
      terra::crs(r) <- "epsg:7589"
      terra::values(r) <- r_vals

      bin <- seg <- equalarea_seg <- rast(r)

      terra::values(seg) <- seg_vals
      names(seg) <- "Sky segments of equal area"
      terra::values(equalarea_seg) <- equalarea_seg_vals
      names(equalarea_seg) <- "Sky segments of equal area"
      terra::values(bin) <- bin_vals
      bin <- binarize_with_thr(bin, 0.5)

      .single(r, bin, seg, equalarea_seg, dist_to_black)
    }
    metrics <- .with_cluster(cores, {
      # Only to avoid note from check, code is OK without this line.
      i <- NA

      foreach(i = 1:nrow(params)) %dopar% {
        .multi(row_vals, col_vals, r_vals,
               bin_vals_list[[params[i, 1]]],
               seg_vals_list[[as.character(params[i, 2])]],
               equalarea_seg_vals,
               params[i, 3])
      }
    })
  } else {
    metrics <- lapply(1:nrow(params),
                       function(i) {
                         .single(r,
                                 bin_list[[params[i, 1]]],
                                 seg_list[[as.character(params[i, 2])]],
                                 equalarea_seg,
                                 params[i, 3])
                         }
                       )
  }
  metrics <- matrix(unlist(metrics), ncol = 2, byrow = TRUE)
  acc_vec <- metrics[,1]
  cov_vec <- metrics[,2]

  acc_n <- normalize_minmax(acc_vec)
  cov_n <- normalize_minmax(cov_vec)

  cam <- (1 - cov_n) * w + acc_n * w
  best <- which.min(cam)

  values <- params[best, ]

  # -------------------------------------------------------------------------
  # Write log and sidecar file if requested
  # -------------------------------------------------------------------------
  if (!is.null(write_log_in)) {

    t1 <- Sys.time()
    timestamp_start <- format(t0, "%Y-%m-%d %H:%M:%S")
    timestamp_end   <- format(t1, "%Y-%m-%d %H:%M:%S")
    elapsed_minutes <- as.numeric(difftime(t1, t0, units = "mins"))

    # bin_list names (fallback si están vacíos)
    bin_names <- names(bin_list)
    if (is.null(bin_names) || all(bin_names == "")) {
      bin_names <- paste0("bin_", seq_along(bin_list))
    }

    # número de celdas del equalarea_seg (indicado por vos)
    equalarea_ncells <- length(unique(terra::values(equalarea_seg)))

    # info del sistema
    sys_i <- Sys.info()
    r_v   <- R.version.string

    # construcción del bloque bin_list_names multilinea y numerado
    bin_names_block <- paste0(
      "  ", seq_along(bin_names), ": ", bin_names,
      collapse = "\n"
    )

    # Log principal ----------------------------------------------------------
    log_lines <- c(
      "Sky-sampling tuning log",
      paste0("Timestamp (start): ", timestamp_start),
      paste0("Timestamp (end):   ", timestamp_end),
      paste0("Elapsed minutes: ", round(elapsed_minutes, 2)),
      "",
      "Test configuration:",
      paste0("  n_cells (equalarea_seg): ", equalarea_ncells),
      paste0("  w: ", w),
      paste0("  logical parallel: ", logical),
      paste0("  parallel: ", parallel),
      paste0("  cores: ", cores),
      "",
      "System information:",
      paste0("  R version: ", r_v),
      paste0("  System: ", sys_i["sysname"], " / ", sys_i["release"]),
      paste0("  Machine: ", sys_i["machine"]),
      "",
      "Parameter grids tested:",
      paste0("  bin_list size: ", length(bin_list)),
      "  bin_list names:",
      bin_names_block,
      paste0("  n_cells_seq: ", paste(n_cells_seq, collapse = ", ")),
      paste0("  dist_to_black_seq: ", paste(dist_to_black_seq, collapse = ", ")),
      "",
      "Selected combination:",
      paste0("  bin_list_index: ", values$bin_list_index),
      paste0("  n_cells: ", values$n_cells),
      paste0("  dist_to_black: ", values$dist_to_black),
      "",
      "Performance metrics:",
      paste0("  accuracy: ", acc_vec[best]),
      paste0("  coverage: ", cov_vec[best]),
      paste0("  cam: ", cam[best]),
      "",
      "Sidecar file containing the full evaluation grid was generated.",
      ""
    )

    writeLines(log_lines, con = paste0(write_log_in, ".txt"))

    # -----------------------------------------------------------------------
    # Sidecar file (CSV compatible con soft de ofimática)
    # -----------------------------------------------------------------------
    df_out <- data.frame(
      params,
      accuracy = acc_vec,
      coverage = cov_vec,
      cam = cam
    )

    sidecar_file <- paste0(write_log_in, "_full_metrics.csv")
    utils::write_csv2(df_out, file = sidecar_file)
  }


  list(n_cells = values$n_cells,
       dist_to_black = values$dist_to_black,
       bin_list_index = values$bin_list_index,
       accuracy = acc_vec[best],
       coverage = cov_vec[best],
       cam = cam[best])
}
