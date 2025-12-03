#' Tune sky-sampling parameters
#'
#' Tune the parameters of [extract_sky_pixels()].
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
#'
#' @inheritParams compute_canopy_openness
#' @inheritParams assess_sampling_uniformity
#' @inheritParams apply_by_direction
#'
#' @details
#' Evaluate combinations of three discrete parameters for sampling sky digital
#' numbers (\eqn{\delta_\text{sky}}) using [extract_sky_points()]:
#' \itemize{
#'   \item the binarized image selected from `bin_list`, provided to
#'     [extract_sky_points()] as argument `bin`;
#'   \item the parameter `n_cells` used to generate a segmentation with
#'     [equalarea_segmentation()], provided to [extract_sky_points()] via `seg`;
#'   \item the minimum allowed distance to black pixels in the binarized image,
#'     provided as argument `dist_to_black` of [extract_sky_points()].
#' }
#'
#' The optimal combination minimizes the Coverage–Accuracy Metric (CAM):
#'
#' \deqn{\mathrm{CAM} = (1 - \mathrm{coverage}) \cdot w \;+\; \mathrm{accuracy} \cdot (1 - w)}
#'
#' CAM balances hemispherical surface coverage and the local accuracy of sampled
#' \eqn{\delta_\text{sky}}, with \eqn{w} acting as a weight that allows the user
#' to prioritize one component over the other.
#'
#' True sky regions typically exhibit smooth radiance gradients, in contrast with
#' bright canopy elements or mixed pixels, which introduce sharper local
#' fluctuations. To capture this contrast, accuracy is estimated using a
#' robust local coefficient of variation computed per sample window as
#' \eqn{\mathrm{CV}_{\text{local}} = \mathrm{MAD}_{3\times3} / \mathrm{median}_{3\times3}},
#' where MAD and median are calculated on a \eqn{3 \times 3} neighborhood.
#' The accuracy reported by the function is the median of these local CV values
#' across the sampled sky points.
#'
#' Coverage is computed with [assess_sampling_uniformity()], i.e., the fraction
#' of cells in `equalarea_seg` containing at least one sky point.
#'
#' Lower CAM values indicate lower local variability and higher angular coverage.
#'
#' @return List with the following components:
#' \describe{
#'   \item{`n_cells`}{Optimal value within `n_cells_seq`.}
#'   \item{`dist_to_black`}{Optimal value within `dist_to_black_seq`.}
#'   \item{`bin_list_index`}{Index of the optimal binarization in `bin_list`.}
#'   \item{`accuracy`}{Numeric. Median of local CV values (MAD / median) at the
#'     sampled sky points. Lower is better.}
#'   \item{`coverage`}{Numeric in \eqn{[0,1]}. Fraction of reference cells with
#'     at least one sampled sky point. Higher is better.}
#'   \item{`cam`}{Numeric. Minimum Coverage–Accuracy Metric (CAM) for the
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
#' sky_points <- extract_sky_points(r,
#'                                  bin_list[[params$bin_list_index]],
#'                                  equalarea_segmentation(z, a, params$n_cells),
#'                                  params$dist_to_black)
#' display_caim(caim$Blue, sky_points = sky_points)
#' }
tune_sky_sampling <- function(r, z, a,
                              equalarea_seg,
                              bin_list,
                              n_cells_seq,
                              dist_to_black_seq,
                              w = 0.5,
                              parallel = TRUE,
                              cores = NULL,
                              logical = TRUE,
                              leave_free = 1) {

  .check_r_z_a_m(r, z, a, m = NULL)
  .assert_single_layer(equalarea_seg)
  .check_vector(n_cells_seq, "integerish", sign = "positive")
  .check_vector(dist_to_black_seq, "integerish", sign = "positive")
  .check_vector(w, "numeric", 1, sign = "nonnegative")
  if (w > 1) (stop("´w´ cannot be greater than one."))
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }

  .single <- function(r, bin, seg, equalarea_seg, dist_to_black) {
    .extract_cv <- function(r, sky_points) {
      #based on extract_d
      cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
      xy <-  terra::xyFromCell(r, cells)

      .fn <- function(x) mad(x) / median(x)
      Map(
        function(x, y) {
          ma <- expand.grid(c(-1,0,1) + x, c(-1,0,1) + y)
          .fn(terra::extract(r, ma, method = "simple", na.rm = TRUE)[,-1])
        },
        xy[,1],
        xy[,2]) %>% as.data.frame() %>% t %>% unname()
    }
    sky_points <- extract_sky_points(r, bin, seg,
                                     dist_to_black = dist_to_black,
                                     method = "per_cell")
    acc <- median(.extract_cv(r, sky_points))
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
  list(n_cells = values$n_cells,
       dist_to_black = values$dist_to_black,
       bin_list_index = values$bin_list_index,
       accuracy = acc_vec[best],
       coverage = cov_vec[best],
       cam = cam[best])
}
