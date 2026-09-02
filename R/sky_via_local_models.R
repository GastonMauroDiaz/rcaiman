#' Sky via local models
#'
#' @description
#' Reconstruct a sky vault by fitting local cone-shaped models to sampled image
#' values.
#'
#' @details
#' This method extends the cone-shaped model proposed by
#' \insertCite{Diaz2018;textual}{rcaiman} by applying it locally across
#' hemispherical space. Neighboring observations are selected around multiple
#' query directions, local cone-shaped models are fitted independently, and
#' their predictions are interpolated to reconstruct a continuous sky vault.
#'
#' Unlike [fit_cone_shaped_model()], which fits a single cone-shaped model to
#' all observations, this approach combines multiple local models, allowing
#' complex sky patterns to emerge without assuming a globally cone-shaped
#' surface.
#'
#' @section Structural redundancy and sampling support:
#'
#' Reconstruction flexibility emerges from the interaction between the density
#' of the query grid (`spacing_deg`) and the spatial distribution of `sky_points`.
#'
#' To assess whether the sampling design allows for true local flexibility, a
#' baseline structural redundancy score is internally computed. This index
#' estimates the average times a single sampling point is reused across
#' neighboring local models under the highest possible local setting.
#'
#' Under ideal, perfectly uniform sampling conditions (\eqn{\mathrm{KL} = 0}),
#' potential structural redundancy (\eqn{R_{\text{ideal}}}) depends only on the
#' ratio between the minimum required points per neighborhood and the available
#' sampling pool:
#'
#' \deqn{R_{\text{ideal}} = \frac{20 \cdot N_{\text{query}}}{N_{\text{points}}}}
#'
#' where \eqn{N_{\text{query}}} is the number of query directions and
#' \eqn{N_{\text{points}}} is the total number of sampling points.
#'
#' In practice, clustered or sparse sampling patterns alter this support. The
#' effective recycling factor is adjusted using the Kullback–Leibler divergence
#' (\eqn{\mathrm{KL}}) obtained from neighborhood support sizes:
#'
#' \deqn{R_{\text{effective}} = R_{\text{ideal}} \cdot (1 + \mathrm{KL})}
#'
#' where \eqn{\mathrm{KL}} is the Kullback–Leibler divergence defined in
#' [assess_sampling_uniformity()], computed here from neighborhood support
#' rather than from occupancy of an equal-area hemispherical tessellation.
#'
#' This quantity is used only as an internal heuristic and should not be
#' interpreted as a formal property of the reconstructed sky surface.
#'
#' A warning is issued when \eqn{R_{\text{effective}}} suggests excessive
#' overlap (values greater than 5), indicating that the data distribution
#' restricts local independence, forcing the reconstruction to behave as if a
#' higher `conicity` was requested.
#'
#' @param col_id character or integer vector of length one. Name or position of
#'   the column in `sky_points` containing the values sampled from the image.
#' @param spacing_deg numeric vector of length one. Angular spacing (degrees)
#'   between query directions where local models are fitted. Smaller values
#'   increase the number of fitted models and have impact on defining what local
#'   means in this context, see *Structural redundancy and sampling support*.
#' @param spacing_px numeric vector of length one. Pixel spacing used to
#'   densify local predictions before reconstructing the final raster. Smaller
#'   values produce denser intermediate predictions.
#' @param conicity numeric vector of length one ranging from 0 to 1. Controls
#'   the tendency of the reconstructed sky surface toward a global cone-shaped
#'   structure. Smaller values fit models using fewer neighboring points and can
#'   capture more complex sky radiance patterns, whereas larger values produce
#'   smoother reconstructions that increasingly resemble a single cone-shaped
#'   surface. The effective flexibility of fitted models also depends on the
#'   spatial distribution of `sky_points`. Highly clustered sampling patterns
#'   may reduce the ability of low `conicity` values to produce highly local
#'   reconstructions.
#'
#' @inheritParams add_directions
#'
#' @return Numeric [terra::SpatRaster-class] with one layer and the same
#' geometry as `z`.
#'
#' @note `spacing_deg` operates as a reconstruction hyperparameter.
#' Smaller values do not necessarily increase reconstruction flexibility because
#' neighboring local models can become highly redundant when sampling support is
#' limited. The redundancy warning can help identify parameter combinations that
#' produce excessive overlap among fitted models.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_jpeg.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry, "exif_jpeg", "delta_t")
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, zenith_col_row = spec$zenith_col_row ,
#'                    horizon_radius = spec$horizon_radius)
#' z <- zenith_image(ncol(caim), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "wagner2001_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`0`)
#' r <- apply_flat_field(caim$Blue, ff)
#' m <- !is.na(z)
#'
#' r <- invert_gamma_correction(r)
#'
#' thr <- apply_by_direction(r, z, a, m)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' seg <- equal_area_segmentation(z, a, 1000)
#' sky_points <- sample_sky_points(r, bin, seg, dist_to_black = 1)
#' sky_points <- add_directions(sky_points, z, a)
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- rem_outliers(
#'    sky_points,
#'    col_id = "Blue",
#'    k = 20,
#'    search_radius = 20,
#'    laxity = 2,
#'    cutoff_side = "left"
#' )
#' cv <- cv_image(r)
#' sky_points <- extract_dn(cv, sky_points, use_window = FALSE)
#' i <- is_outlier(sky_points$cv, sky_points$cv, laxity = 1, cutoff_side = "right")
#' sky_points <- sky_points[!i,]
#'
#' plot(m)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky <- sky_via_local_models(z, a, sky_points, col_id = "Blue",
#'                              spacing_deg = 25, conicity = 0.5)
#' plot(sky)
#' plot(r/sky)
#'
#' thr <- thr_mblt(sky, -7.8, 0.95 * 0.5)
#' plot(binarize_with_thr(r, thr))
#' }
sky_via_local_models <- function(z,
                                 a,
                                 sky_points,
                                 col_id,
                                 conicity = 0,
                                 spacing_deg = 10,
                                 spacing_px = 5
) {
  col_id <- .resolve_required_column(col_id, sky_points)
  .check_columns(sky_points,
                 required_columns = c("row", "col", "zen", "azi", col_id),
                 reserved_columns = NULL)
  .check_vector(conicity, "numeric", 1, sign = "nonnegative")

  if (conicity > 1) {
    stop("`conicity` cannot be greater than one.")
  }

  k_cs <- round(
    20 + conicity * (nrow(sky_points) - 20)
  )

  query_indices <- fibonacci_indices(z, a, spacing_deg)

  query_indices$seed <- seq_len(nrow(query_indices))
  nn_density <- interpolate_spherical(
    query_indices,
    sky_points,
    col_id = "seed",
    k = 1,
    search_radius = 180,
    rule = "any"
  )
  counts <- tapply(nn_density$seed, nn_density$seed, length)

  p <- counts / sum(counts)

  N <- length(counts)

  kl <- sum(
    p[p > 0] * log(p[p > 0] * N)
  )

  # Medimos la redundancia estructural pura de los datos (peor escenario local k = 20)
  baseline_redundancy <- 20 * nrow(query_indices) / nrow(sky_points)
  recycling_factor <- baseline_redundancy * (1 + kl)

  # Un umbral > 5 indica que cada punto se recicla excesivamente, destruyendo la independencia local
  if (recycling_factor > 5) {
    warning(
      paste0(
        "High data redundancy detected: each sampling point is recycled approx. ",
        round(recycling_factor, 1), " times across local models. ",
        "This reduces local flexibility, forcing the reconstruction to behave ",
        "as if a higher `conicity` was requested. ",
        "Consider increasing `spacing_deg` to restore local independence."
      )
    )
  }

  nn <- nn_spherical(sky_points, query_points,
                     k = k_cs,
                     search_radius = NULL,
                     rule = "any")
  models <- apply_locally(
    sky_points,
    query_points,
    nn,
    fun = function(df, ...) {
      fit_cone_shaped_model(df, col_id = col_id)
    },
    parallel = FALSE
  )
  pred <- lapply(seq_along(models), function(i) {
    if (length(models[[i]]) == 1) return(NA)
    model <- models[[i]]
    model$fun(query_points[i, "zen"], query_points[i, "azi"], units = "radian")
  }) %>% unlist()
  eval_points <- cbind(query_points, pred)
  eval_points <- eval_points[!is.na(pred),]
  raster_grid_pts <- raster_grid_in_circle(z, a, spacing_px)
  pts <- interpolate_spherical(
    eval_points,
    raster_grid_pts,
    col_id = "pred",
    k = 3,
    p = 2,
    search_radius = spacing_deg * 2,
    rule = "any"
  )
  sky <- interpolate_planar(
    pts,
    col_id = "pred",
    z,
    k = 3,
    p = 2,
    search_radius = spacing_px * 2
  )
  sky[is.na(z)] <- 0
  sky
}
