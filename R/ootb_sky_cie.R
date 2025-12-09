#' Out-of-the-box CIE sky model and raster
#'
#' @description
#' Fit and validate a CIE general sky model to radiance sampled on canopy
#' gaps and return the predicted raster.
#'
#' @note
#' This function is part of a paper under preparation.
#'
#' @param twilight logical vector of length one. When `TRUE`, test sun zenith
#'   angles from 90 to 96 deg (civil twilight), plus five additional zenith
#'   angles equally spaced between the estimated sun zenith angle and the
#'   horizon. This is required because [estimate_sun_angles()] with
#'   `method = "assume_obscured"` may incorrectly identify the center of the
#'   visible circumsolar region as the solar disk.
#'
#' @inheritParams tune_sky_sampling
#' @inheritParams compute_canopy_openness
#' @inheritParams fit_cie_model
#' @inheritParams apply_by_direction
#'
#' @return List with:
#' \describe{
#'   \item{`rr_raster`}{numeric [terra::SpatRaster-class]. Predicted relative radiance.}
#'   \item{`model`}{list returned by [fit_cie_model()]. The optimal fit.}
#'   \item{`params`}{list returned by [tune_sky_sampling()].}
#'   \item{`sky_points`}{`data.frame` with columns `row` and `col`. Locations of
#'     sky points.}
#'   \item{`sun_row_col`}{`data.frame` with the estimated sun‑disk position in
#'     image coordinates.}
#'   \item{`optimal_start`}{`custom_sky_coef` with the estimated sun‑disk position in
#'     image coordinates.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' com <- complementary_gradients(caim)
#' thr <- thr_twocorner(com$red_cyan[])
#' bin <- binarize_with_thr(com$red_cyan, thr$uc)
#' bin <- rem_isolated_black_pixels(bin)
#' bin_list <- c(bin_list, bin)
#' bin_list <- c(bin_list, rem_small_gaps(bin))
#'
#' set.seed(7)
#' sky_cie <- ootb_sky_cie(r, z, a, m,
#'                         bin_list,
#'                         n_cells_seq,
#'                         dist_to_black_seq,
#'                         method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
#'                         twilight = FALSE,
#'                         std_sky_no = 12,
#'                         optim_zenith_dn = TRUE,
#'                         parallel = TRUE)
#'
#' sky_cie$rr_raster
#' plot(sky_cie$rr_raster)
#' plot(sky_cie$model$rr$sky_points$pred, sky_cie$model$rr$sky_points$rr,
#'      xlab = "Predicted relative radiance", ylab = "Observed relative radiance")
#' abline(0,1)
#'
#' ratio <- r/sky_cie$rr_raster/sky_cie$model$rr$zenith_dn
#' plot(ratio)
#' plot(select_sky_region(ratio, 0.95, 1.05))
#' plot(select_sky_region(ratio, 1.05, 100))
#'
#' display_caim(caim, sky_points = sky_cie$sky_points, sun_row_col = sun)
#' }
ootb_sky_cie <- function(r, z, a, m,
                         bin_list,
                         n_cells_seq,
                         dist_to_black_seq,
                         w = 0.5,
                         method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
                         twilight = FALSE,
                         custom_sky_coef = NULL,
                         std_sky_no = NULL,
                         general_sky_type = NULL ,
                         optim_zenith_dn = FALSE,
                         parallel = TRUE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1
                         ) {

  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }
  #more basic checks are handled by the functions called below


# Tune --------------------------------------------------------------------

  params <- tune_sky_sampling(
    r,
    z,
    a,
    equalarea_segmentation(z, a, 100),
    bin_list = bin_list,
    n_cells_seq = n_cells_seq,
    dist_to_black_seq = dist_to_black_seq,
    w = w,
    parallel = parallel,
    cores = cores,
    logical = logical,
    leave_free = leave_free
  )

  seg <- equalarea_segmentation(z, a, params$n_cells)
  bin <- bin_list[[params$bin_list_index]]
  dist_to_black <- params$dist_to_black

# Get marks ---------------------------------------------------------------

  sun_angles <- estimate_sun_angles(r, z, a, bin, seg, angular_radius_sun = 30,
                                    method = "assume_obscured")

  sun_angles2 <- estimate_sun_angles(r, z, a, bin, NULL, NULL,
                                     method = "assume_veiled")

  sky_points <- sample_sky_points(r, bin, seg, dist_to_black)


# Apply filters to sky points ---------------------------------------------
  sky_points <- rem_nearby_points(sky_points, r, min_dist = 3, space = "planar")

  i <- apply_local_spherical(
    cbind(sky_points, extract_cv(r, sky_points)),
    NULL,
    z,
    a,
    k = 20,
    angular_radius = 45,
    rule = "all",
    fun = function(df) {
      !is_outlier(df[1, 3], df[-1, 3])
    },
    parallel = FALSE
  )
  #query points with nn < k
  i$output[is.na(i$output)] <- TRUE

  sky_points <- rem_outliers(
    sky_points[i$output, ],
    r,
    z,
    a,
    k = 20,
    angular_radius = 45,
    laxity = 2,
    cutoff_side = "right",
    use_window = TRUE,
    detrend = FALSE
  )

  sky_points <- rem_outliers(
    sky_points,
    r,
    z,
    a,
    k = 20,
    angular_radius = 30,
    laxity = 3,
    cutoff_side = "left",
    use_window = TRUE,
    detrend = TRUE
  )

  equalarea_seg <- equalarea_segmentation(z, a, 200)
  kl <- lapply(1:10, function(i) {
    sky_points <-  rem_nearby_points(sky_points, NULL, z, a, i, space = "spherical")
    assess_sampling_uniformity(sky_points, equalarea_seg)$kl
  })
  sky_points <- rem_nearby_points(sky_points, r, z, a, min_dist = which.min(kl), space = "spherical")

# Fit ---------------------------------------------------------------------

  rr <- extract_rr(r, z, a, sky_points, no_of_points = 3, use_window = TRUE)

  if (twilight) {
    civic_twilight <- c(seq(sun_angles["z"], 90, length = 5), seq(91, 96, 1))
    suns <- lapply(seq_along(civic_twilight),
                   function(i) c(z = civic_twilight[i], sun_angles["a"]))
    suns[[length(suns) + 1]] <- sun_angles2
  } else {
    suns <- list(sun_angles, sun_angles2)
  }

  models <- if (parallel && twilight) {
    # Only to avoid note from check, code is OK without this line.
    sun <- NA

    .with_cluster(cores, {
      foreach(sun = suns) %dopar% {
        fit_cie_model(rr, sun, custom_sky_coef = custom_sky_coef,
                      std_sky_no = NULL,
                      general_sky_type = NULL ,
                      optim_zenith_dn = optim_zenith_dn,
                      method = method)
      }
    })
  } else {
    lapply(suns, function(sun) {
      fit_cie_model(rr, sun, custom_sky_coef = custom_sky_coef,
                    std_sky_no = std_sky_no,
                    general_sky_type = general_sky_type,
                    optim_zenith_dn = optim_zenith_dn,
                    method = method)
    })
  }

  metrics <- lapply(seq_along(models),
                    function(i) models[[i]]$metric) %>% unlist
  i <- which.min(metrics)

  model <- models[[i]]

# Optim sun ---------------------------------------------------------------

  model <- optim_sun_angles(model, method = method)

# output ------------------------------------------------------------------

  .get_sky_cie <- function(z, a, model) {
    cie_image(z, a, model$sun_angles,  model$coef)
  }

  sun_row_col <- row_col_from_zenith_azimuth(z, a,
                                             model$sun_angles["z"],
                                             model$sun_angles["a"])

  list(rr_raster = .get_sky_cie(z, a, model),
       model = model,
       params = params,
       sky_points = model$rr$sky_points[, c("row", "col")],
       sun_row_col = sun_row_col
       )
}
