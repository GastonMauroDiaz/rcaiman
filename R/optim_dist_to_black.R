#' Optimize minimum distance to black pixels
#'
#' @description
#' Estimate an optimal buffer (`dist_to_black`) to keep sampled sky points away
#' from candidate canopy pixels (black pixels).
#'
#' @details
#' The heuristic seeks the largest buffer that still yields uniform angular
#' coverage. It iteratively decreases `dist_to_black` while monitoring the
#' percentage of 30 deg sky‑grid cells covered by sampled points. If coverage
#' is low, the buffer is relaxed (and may be removed). This balances border
#' avoidance with representativeness across the sky vault.
#'
#' @inheritParams skygrid_segmentation
#' @inheritParams extract_sky_points
#' @inheritParams compute_canopy_openness
#'
#' @return Numeric vector of length one to be passed as `dist_to_black` to
#'   [extract_sky_points()], or `NULL` if no buffer is advised.
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' bin <- binarize_by_region(r, ring_segmentation(z, 15), "thr_isodata") &
#'   select_sky_region(z, 0, 88)
#' g <- skygrid_segmentation(z, a, 10, first_ring_different = TRUE)
#'
#' dist_to_black <- optim_dist_to_black(r, z, a, m, bin, g)
#' dist_to_black
#'
#' bin <- grow_black(bin, 11)
#' plot(bin)
#' dist_to_black <- optim_dist_to_black(r, z, a, m, bin, g)
#' dist_to_black
#'
#' bin <- grow_black(bin, 21)
#' plot(bin)
#' dist_to_black <- optim_dist_to_black(r, z, a, m, bin, g)
#' dist_to_black
#' }
#' @export
optim_dist_to_black <- function(r, z, a, m, bin, g,
                                parallel = FALSE,
                                cores = NULL,
                                logical = TRUE,
                                leave_free = 1) {
  .check_r_z_a_m(r, z, a, m, r_type = "single")
  .assert_logical_mask(bin)
  .assert_single_layer(g)
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }

  sky_points <- extract_sky_points(r, bin, g, dist_to_black = NULL)
  equalarea_seg <- equalarea_segmentation(z, a, 32)
  kl <- assess_sampling_uniformity(sky_points, equalarea_seg)

  # .get_kl <- function(i) {
  #
  # }

  # if (parallel) {
  #   acc <- .with_cluster(cores, {
  #     i_chunks <- split(seq_len(n_directions), 1:cores) %>% suppressWarnings()
  #
  #     # Only to avoid note from check, code is OK without this line.
  #     j <- NA
  #
  #     foreach::foreach(j = 1:cores,
  #                      .combine = rbind,
  #                      .packages = c("terra")) %dopar% {
  #                        do.call(rbind, lapply(i_chunks[[j]],
  #                                              .extract_fov_n_compute_fun))
  #                      }
  #   })
  # } else {
  #   acc <- do.call(rbind, lapply(seq_len(n_directions),
  #                                .extract_fov_n_compute_fun))
  # }


  d <- seq(1, 11, 2) %>% rev()
  kls <- c()
  for (i in seq_along(d)) {
    sky_points <- extract_sky_points(r, bin, g,
                                     dist_to_black = d[i])
    kls <- c(kls, assess_sampling_uniformity(sky_points, equalarea_seg))
  }
browser()



  g30[!m] <- 0
  dist_to_black <- 11
  sampling_pct <- 0
  while (sampling_pct < 100 & dist_to_black > 3) {
    dist_to_black <- dist_to_black - 2
    sky_points <- extract_sky_points(r, bin, g,
                                     dist_to_black = dist_to_black)
    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
      (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
  }
  if (sampling_pct < 75) {
    dist_to_black <- 1
    sky_points <- extract_sky_points(r, bin, g,
                                     dist_to_black = dist_to_black)
    v <- cellFromRowCol(r, sky_points$row, sky_points$col) %>%
      xyFromCell(r, .) %>% vect()
    sampling_pct <- (extract(g30, v)[,2] %>% unique() %>% length()) /
      (unique(g30)[,1] %>% length() %>% subtract(1)) * 100
  }
  if (sampling_pct < 50) {
    dist_to_black <- NULL
  }
  dist_to_black
}
