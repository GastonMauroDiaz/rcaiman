#' Sky via local models
#'
#' @description
#' Reconstruct a sky image by fitting local cone-shaped models to sampled
#' values and interpolating their predictions.
#'
#' @details
#' This method extends the cone-shaped model proposed by
#' \insertCite{Diaz2018;textual}{rcaiman} by applying it locally across the sky.
#' Unlike [fit_cone_shaped_model()], which fits one model to all observations,
#' it allows different regions to follow different cone-shaped models.
#' `conicity` controls the transition from local fits to a global fit without
#' requiring the user to choose a number of observations for each model.
#'
#' All calculations of angular neighborhoods and cone-shaped models use
#' external directions, relative to zenith and north. Supply the acquisition
#' orientation to [direction_image()]. Raster angles are in degrees and table
#' angles are in radians, following [add_directions()].
#'
#' Query directions are generated on a Fibonacci lattice over the upper
#' hemisphere with approximate separation `spacing_deg`. They remain continuous
#' and are generated independently of image coverage and `m`, including outside
#' the observed field of view. Each query uses its nearest observations without
#' an angular search limit. Select suitable sky observations before calling this
#' function; the output mask does not select observations for fitting.
#'
#' Each local fit uses the `"zenith_n_azimuth"` method of
#' [fit_cone_shaped_model()] and is evaluated only at its query direction.
#' Failed fits return no prediction, with the warning from that function.
#' Predictions are interpolated in spherical space onto selected cells of an
#' intermediate image grid spaced by `spacing_px`. These values are then
#' interpolated in planar space onto the final grid. Both steps use inverse
#' distance weighting with up to three neighbors and power two. Search radii
#' are `min(180, 2 * spacing_deg)` degrees and `2 * spacing_px` pixels,
#' respectively. Neither step uses background filling or clips predicted values.
#'
#' The final planar step is an approximation that depends on image projection
#' and resolution. Decreasing `spacing_px` makes the intermediate grid denser;
#' compare reconstructions to assess its effect for the image in question.
#'
#' @section Neighborhoods and sampling support:
#' Each fit uses the nearest observations in angular distance, without a fixed
#' angular radius. For \eqn{N} observations, the requested number of neighbors
#' is `round(20 + conicity * (N - 20))`. With at least 20 observations,
#' `conicity = 0` uses 20 neighbors and `conicity = 1` uses all observations.
#'
#' The selected observations lie in a spherical cap centered on the query
#' (see [nn_spherical()]), whose angular radius is the distance to the farthest
#' selected neighbor. Here that radius results from selecting a number of
#' observations, rather than being specified in advance.
#'
#' The angular extent of each neighborhood depends on the distribution of
#' observations: densely sampled regions can supply the requested neighbors
#' nearby, whereas sparsely sampled regions require more distant observations.
#' A low `conicity` therefore does not guarantee a small angular neighborhood.
#'
#' Reducing `spacing_deg` increases the number of query directions, but does
#' not add observations or guarantee finer reconstructed detail. Nearby queries
#' can use overlapping or identical sets of observations. Identical sets give
#' the same fitted model, evaluated at different directions. Such overlap is
#' expected and does not by itself indicate a problem with the reconstruction.
#'
#' With `conicity = 1`, every query uses the same global fit. The output still
#' passes through the two interpolation steps, so it is not generally identical
#' to evaluating that model at every pixel. `spacing_px` controls the density
#' of the intermediate image grid; it does not change the observations used
#' for fitting.
#'
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, in that order. Defines the output grid
#'   and the external direction of each pixel.
#' @param sampling_directions data frame with numeric columns `zenith` and
#'   `azimuth` in radians and the observed values selected by `col_id`.
#'   Required columns must be finite; rows are not removed automatically.
#'   At least 20 observations are needed. Additional columns, including cell
#'   indices, are ignored. See [fit_cone_shaped_model()] for angle ranges and
#'   sample requirements.
#' @param col_id character or integer vector of length one. Name or position
#'   of the column containing observed values in `sampling_directions`.
#' @param conicity numeric vector of length one ranging from 0 to 1. Controls
#'   the tendency toward a global cone-shaped model by increasing the number
#'   of observations used in each fit. Zero uses the 20 nearest observations;
#'   one uses all observations. See *Neighborhoods and sampling support*.
#' @param spacing_deg positive numeric vector of length one, at most 180.
#'   Approximate separation in degrees between query directions where local
#'   models are fitted and evaluated. Also sets the spherical interpolation
#'   radius as described in Details.
#' @param spacing_px positive integer vector of length one. Spacing of rows
#'   and columns in the intermediate image grid, starting at row and column
#'   one. Also sets the planar interpolation radius as described in Details.
#' @param m optional logical, NA-free [terra::SpatRaster-class] on the same
#'   grid as `d`. Limits the output and intermediate grid to selected cells.
#'   With `NULL`, all cells with both angles available and zenith angle at
#'   most 90 degrees are selected. These angular restrictions also apply when
#'   a mask is supplied. An undefined azimuth, including at zenith, is excluded.
#'
#' @return Numeric [terra::SpatRaster-class] with one layer named `sky` and the
#'   same grid as `d`. Cells outside the selected sky region or lacking
#'   interpolation support are `NA`. With fewer than 20 observations, returns
#'   an all-`NA` raster with a warning. An empty selected intermediate grid or
#'   failure of all local fits also leaves the result entirely `NA`.
#'
#' @seealso [fit_cone_shaped_model()], [cone_shaped_image()],
#'   [interpolate_spherical()], [interpolate_planar()]
#' @references \insertAllCited{}
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
#' d <- direction_image(g, elevation = 60)
#' m <- !is.na(d$zenith) & !is.na(d$azimuth) & d$zenith <= 90
#'
#' # Illustrative sky values. With measurements, select sky observations from
#' # a radiometrically corrected image before fitting.
#' r <- 1000 + 200 * (d$zenith * pi / 180)^2 +
#'   100 * sin(d$azimuth * pi / 180)^2
#' names(r) <- "Blue"
#' sampling_indices <- expand.grid(col = seq(1, 80, 4), row = seq(1, 60, 4))
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' samples <- extract_dn(r, sampling_indices, use_window = FALSE)
#' samples <- add_directions(samples, d = d)
#' sky <- sky_via_local_models(d, samples, "Blue", conicity = 0.3,
#'                             spacing_deg = 15, spacing_px = 3, m = m)
#' plot(sky)
#' global <- sky_via_local_models(d, samples, "Blue", conicity = 1,
#'                                spacing_deg = 15, spacing_px = 3, m = m)
#' plot(sky - global)
sky_via_local_models <- function(d,
                                 sampling_directions,
                                 col_id,
                                 conicity = 0,
                                 spacing_deg = 10,
                                 spacing_px = 5,
                                 m = NULL) {
  z <- .sky_segmentation_layer(d, "zenith")
  a <- .sky_segmentation_layer(d, "azimuth")
  .check_vector(conicity, "numeric", 1, sign = "nonnegative")
  .check_vector(spacing_deg, "numeric", 1, sign = "positive")
  .check_vector(spacing_px, "integerish", 1, sign = "positive")
  if (conicity > 1) stop("`conicity` cannot be greater than one.", call. = FALSE)
  if (spacing_deg > 180) {
    stop("`spacing_deg` cannot exceed 180 degrees.", call. = FALSE)
  }
  .check_columns(sampling_directions, c("zenith", "azimuth"))
  if (is.numeric(col_id)) .check_vector(col_id, "integerish", 1, sign = "positive")
  col_id <- .resolve_required_column(col_id, sampling_directions)
  for (column in unique(c("zenith", "azimuth", col_id))) {
    .check_vector(sampling_directions[[column]], "numeric", name = column)
  }
  .spherical_neighbor_vectors(sampling_directions, c("zenith", "azimuth"))
  selected <- !is.na(z) & !is.na(a) & z <= 90
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
    selected <- selected & m
  }
  sky <- terra::rast(z)
  terra::values(sky) <- NA_real_
  names(sky) <- "sky"
  if (nrow(sampling_directions) < 20) {
    warning("At least 20 observations are required to fit the model.", call. = FALSE)
    return(sky)
  }
  grid <- expand.grid(col = seq(1, terra::ncol(z), by = spacing_px),
                      row = seq(1, terra::nrow(z), by = spacing_px))
  grid <- rem_unselected(grid, selected)
  if (!nrow(grid)) return(sky)
  grid <- add_directions(grid, d = d)

  # This lattice describes the sky, independently of the instrument frame.
  queries <- as.data.frame(.cartesian_to_directions(.fibonacci(spacing_deg)))
  names(queries) <- c("zenith", "azimuth")
  k <- round(20 + conicity * (nrow(sampling_directions) - 20))
  neighbors <- nn_spherical(sampling_directions, queries, k = k,
                             search_radius = NULL, parallel = FALSE)
  predictions <- apply_locally(sampling_directions, queries, neighbors,
    fun = function(data, distance, i) {
      fit <- fit_cone_shaped_model(data, col_id = col_id, return = "all")
      if (is.null(fit)) return(NA_real_)
      as.numeric(stats::predict(fit$model, newdata = queries[i, , drop = FALSE]))
    })
  names(predictions)[names(predictions) == "output"] <- "prediction"
  if (all(is.na(predictions$prediction))) return(sky)
  grid <- interpolate_spherical(predictions, grid, col_id = "prediction",
    k = 3, p = 2, search_radius = min(180, spacing_deg * 2), rule = "any")
  sky <- interpolate_planar(add_coord(grid), col_id = "prediction", r = z,
    k = 3, p = 2, search_radius = spacing_px * 2, parallel = FALSE)
  sky[!selected] <- NA_real_
  names(sky) <- "sky"
  sky
}
