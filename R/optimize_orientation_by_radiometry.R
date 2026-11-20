#' Optimize acquisition orientation by radiometric agreement
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Refine an acquisition orientation by comparing its sky samples with an
#' unobstructed reference image. The radiometric scale and reference orientation
#' remain fixed while the sampled viewing directions are rotated.
#'
#' @param sampling_directions data frame with numeric `off_axis` and
#'   `internal_azimuth` columns in radians, and measured values in `col_id`.
#'   These are sky samples from the rover acquisition (see *Acquisition and
#'   preparation*).
#'   Obtain the directions with [add_directions()] or
#'   [internal_directions_from_coord()]. All rows must have complete directions
#'   and finite measurements. Supply at least three rows. Azimuth may be `NA`
#'   in the optical-axis direction or its opposite direction.
#' @param col_id character vector of length one naming the measurement column.
#' @param reference single-layer numeric [terra::SpatRaster-class] containing
#'   the sky reference prepared from the base acquisition. Its rows and columns
#'   must correspond to
#'   `reference_geometry`. Image coordinates follow [y_down_from_y()]; raster
#'   extent and CRS are not used. `NA` cells are unavailable for interpolation.
#' @param reference_geometry radial geometry list following
#'   [set_internal_geometry()]. The source samples can come from either
#'   geometric model, but the reference must use `"radial_projection"`.
#' @param initial_orientation,reference_orientation named lists with
#'   `sensor_rotation`, `azimuth_rotation`, and `elevation`, in degrees, using
#'   the convention of [direction_image()]. `initial_orientation` is the
#'   starting orientation of the rover acquisition; `reference_orientation`
#'   describes the reference and is held fixed. Both default to an upward
#'   acquisition with sensor top toward north.
#' @param reference_scale positive numeric vector of length one. Multiplier
#'   converting the base reference values to the radiometric scale of the rover
#'   measurements. Use 1 when the prepared data already share a common scale,
#'   including after applying cross-calibration. See *Acquisition and preparation*
#'   for independently calibrated instruments.
#' @param max_rotation positive numeric vector of length one, at most 180.
#'   Maximum rotation from `initial_orientation`, in degrees. This bounds the
#'   angle of the complete three-dimensional rotation, not each orientation
#'   argument separately.
#' @param n_starts positive integer-like vector of length one. Number of
#'   deterministic starting rotations, including the initial orientation.
#' @param method optimization method passed to [stats::optim()]:
#'   `"Nelder-Mead"`, `"BFGS"`, or `"CG"`.
#' @param max_iter positive integer-like vector of length one. Maximum number
#'   of iterations of [stats::optim()] per starting rotation.
#'
#' @section Acquisition and preparation:
#' In a base-rover workflow, the base instrument stays at a fixed location while
#' the rover moves between sampling points under the canopy. The instruments
#' need not be identical. Choose a nearby base location with enough sky
#' visibility to prepare the reference. This can be inside the forest, where the
#' canopy is more open and foliage less dense; then, prepare the base reference
#' with CIE modeling and spherical interpolation where needed to replace canopy
#' obstructions.
#'
#' Coordinate each rover acquisition with a base acquisition so that both
#' describe the same sky radiance distribution. Two operators can coordinate
#' their exposures by a spoken signal, or use remote triggers.
#'
#' The reference surface must represent unobstructed sky in the directions used
#' for comparison although the base photograph itself may contain foliage. This
#' concentrates sky reconstruction at a location chosen for favorable sampling
#' conditions. Sky samples from the rover are then used to refine its
#' orientation against that reference. Supply the prepared surface as
#' `reference`; the function treats measured and reconstructed values alike.
#'
#' Apply geometric calibration and flat-field correction before fitting.
#' Establish the radiometric scale by either of these two routes:
#'
#' 1. **Cross-calibrated instruments.** Apply the calibration obtained with
#'    [rcaiman.registry::crosscalibrate_radiometry()] so that the prepared base
#'    and rover data share a common radiometric scale. Use `reference_scale = 1`.
#' 2. **Independently calibrated instruments.** Calibrate each instrument with
#'    the simple method in [rcaiman.registry::calibrate_radiometry()]. This
#'    corrects each instrument's response relative to its optical axis but does
#'    not equalize the scales of the two instruments. Determine a multiplier
#'    from comparable measurements of a common radiometric reference. For
#'    example, `reference_dn["rover"] / reference_dn["base"]` converts base
#'    values to the rover scale and can be supplied as `reference_scale`.
#'
#' The radiometric corrections and relative exposures must match those used
#' to establish the common scale or multiplier. If the multiplier has already
#' been applied to `reference`, use `reference_scale = 1`. The function holds
#' this scale fixed; it does not estimate a gain or intercept.
#'
#' Neither perfect levelling nor identical azimuthal orientation is required.
#' Supply an approximate rover orientation and an appropriate search range;
#' the declared reference orientation remains fixed.
#'
#' Select sky samples before calling the function. Their distribution and purity
#' affect the orientation estimate. Every supplied row participates throughout
#' the search. The initial orientation must place all samples where the
#' reference can be interpolated. Rows that do not meet this condition are
#' reported by an error. Leave enough coverage around the samples for the
#' intended rotation. Candidates losing support for any sample are rejected.
#'
#' @details
#' The scientific method is also experimental. Its effectiveness in operational
#' rover/base acquisitions is still being evaluated. A better radiometric fit
#' alone does not establish a more accurate acquisition orientation.
#'
#' Each candidate rotation converts the internal sample directions to external
#' directions, then into the fixed reference frame. The reference radial model
#' maps these directions to continuous image coordinates. Bilinear interpolation of
#' reference pixel centers supplies predicted measurements, multiplied by
#' `reference_scale`. All contributing pixels must be finite; interpolation
#' does not bridge missing cells or extend beyond the outer pixel centers.
#' No full image reprojection is needed during fitting.
#'
#' The objective is mean squared deviation between measured and predicted
#' values, with equal weight for each sample. Samples and radiometric scale
#' remain fixed as described in "Acquisition and preparation".
#'
#' The search uses a rotation vector: its direction defines a rotation axis
#' in the external east-north-up frame, and its length defines the angle.
#' This rotation is applied to the initial orientation. It avoids optimizing
#' the three orientation arguments directly near a vertical optical axis,
#' where sensor rotation and horizontal rotation are not separately identifiable.
#'
#' The first start is the initial orientation. The next six are positive and
#' negative rotations around each external axis, at half `max_rotation`.
#' Additional starts use axes distributed by a Fibonacci construction at the
#' same angular distance. Each valid start is refined continuously with
#' [stats::optim()]. The best finite candidate encountered is retained, even
#' when an optimization stops without convergence. Invalid starts are reported
#' and skipped. Equal costs retain the earlier candidate, so a uniform reference
#' leaves the initial orientation unchanged. No random sampling is involved.
#'
#' The result optimizes radiometric agreement; it need not uniquely determine
#' physical orientation. A uniform sky provides no directional information,
#' and an axially symmetric sky cannot resolve rotation about its symmetry
#' axis. Calibration errors and reconstruction of the reference also affect
#' the result. Inspect the fitted comparison and, when useful, reproject the
#' acquisitions to a common destination for visual refinement.
#'
#' @return A plain list containing:
#'   \describe{
#'     \item{`orientation`}{named list of the three arguments for
#'       [direction_image()], in degrees. At a vertical optical axis,
#'       `sensor_rotation` is set to zero, as in [fit_image_orientation()].}
#'     \item{`samples`}{input table in original row order, with fitted external
#'       `zenith` and `azimuth` in radians, `reference_x` and `reference_y_down`
#'       in pixels, `predicted` on the measurement scale, and `residual`
#'       (measurement minus prediction). Columns with these names are replaced.}
#'     \item{`metric`, `initial_metric`}{final and initial mean squared
#'       deviations, in squared measurement units.}
#'     \item{`rotation_angle`}{angle in degrees of the rotation from the initial
#'       orientation to the selected orientation.}
#'     \item{`attempts`}{data frame with each starting rotation vector in degrees
#'       (`east`, `north`, `up`), `initial_metric`, best encountered `metric`,
#'       `convergence` code from [stats::optim()], and `message`. Convergence is
#'       `NA` for an invalid start or an optimization error.}
#'   }
#'
#' @seealso [fit_image_orientation()], [direction_image()], [reproject_image()],
#'   [fit_cie_model()], [interpolate_spherical()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(50, 50),
#'   normalization_radius = 50,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(100, 100), geometry)
#' # A synthetic sky with spatial structure and a known acquisition orientation.
#' sky <- function(zenith, azimuth) {
#'   100 + 20 * cos(zenith) + 15 * sin(zenith)^2 * cos(2 * azimuth) +
#'     10 * sin(zenith)^3 * sin(3 * azimuth)
#' }
#' d <- direction_image(g)
#' reference <- sky(d$zenith * pi/180, d$azimuth * pi/180)
#' samples <- fibonacci_directions(8)
#' samples <- samples[samples$off_axis < 55 * pi/180, ]
#' known <- external_directions_from_internal(samples,
#'   sensor_rotation = 20, azimuth_rotation = 22, elevation = 87)
#' samples$Blue <- sky(known$zenith, known$azimuth)
#' fit <- optimize_orientation_by_radiometry(
#'   samples, "Blue", reference, geometry)
#' fit$orientation
#' c(before = fit$initial_metric, after = fit$metric)
#' oriented_d <- do.call(direction_image, c(list(g = g),
#'                                         fit$orientation))
optimize_orientation_by_radiometry <- function(
    sampling_directions, col_id, reference, reference_geometry,
    initial_orientation = list(sensor_rotation = 0, azimuth_rotation = 0,
                               elevation = 90),
    reference_orientation = list(sensor_rotation = 0, azimuth_rotation = 0,
                                 elevation = 90),
    reference_scale = 1, max_rotation = 10, n_starts = 7,
    method = "Nelder-Mead", max_iter = 200) {
  .check_columns(sampling_directions,
                 required_columns = c("off_axis", "internal_azimuth"))
  .check_vector(col_id, "character", 1)
  .check_columns(sampling_directions, required_columns = col_id)
  if (nrow(sampling_directions) < 3L) {
    stop("Supply at least three sky samples.", call. = FALSE)
  }
  rays <- .orientation_unit_vectors(sampling_directions,
                                    "off_axis", "internal_azimuth")
  measured <- sampling_directions[[col_id]]
  .check_vector(measured, "numeric", name = col_id)
  .assert_single_layer(reference)
  .check_radial_destination(reference_geometry)
  .check_vector(reference_scale, "numeric", 1, sign = "positive")
  .check_vector(max_rotation, "numeric", 1, sign = "positive")
  if (max_rotation > 180) stop("`max_rotation` must be at most 180 degrees.",
                              call. = FALSE)
  .check_vector(n_starts, "integerish", 1, sign = "positive")
  .check_vector(max_iter, "integerish", 1, sign = "positive")
  method <- match.arg(method, c("Nelder-Mead", "BFGS", "CG"))
  initial_frame <- .radiometric_orientation_frame(initial_orientation,
                                                 "initial_orientation")
  reference_frame <- .radiometric_orientation_frame(reference_orientation,
                                                   "reference_orientation")
  pixels <- terra::as.matrix(reference, wide = TRUE)
  center <- reference_geometry[["center_of_vision"]]
  limit <- reference_geometry[["max_off_axis_angle"]]
  if (is.null(limit)) limit <- 90
  max_radians <- max_rotation * pi/180

  evaluate <- function(offset) {
    external <- rays %*% t(.orientation_rotation(offset * max_radians) %*%
                            initial_frame)
    internal <- external %*% reference_frame
    angles <- .angles_from_ray(internal[, 1], internal[, 2], internal[, 3],
                               c("off_axis", "internal_azimuth"))
    radius <- reference_geometry[["normalization_radius"]] *
      .calc_relative_radius(angles[, 1], reference_geometry[["parameters"]])
    x <- center[1] + radius * sin(angles[, 2])
    y <- center[2] - radius * cos(angles[, 2])
    on_axis <- angles[, 1] == 0
    x[on_axis] <- center[1]
    y[on_axis] <- center[2]
    outside <- angles[, 1] > limit * pi/180 + 16 * .Machine$double.eps |
      angles[, 1] == pi
    x[outside] <- y[outside] <- NA_real_
    predicted <- reference_scale * .orientation_bilinear(pixels, x, y)
    list(metric = if (all(is.finite(predicted)))
           mean((measured - predicted)^2) else Inf,
         predicted = predicted, x = x, y = y, external = external)
  }
  initial <- evaluate(c(0, 0, 0))
  if (!is.finite(initial$metric)) {
    bad <- which(!is.finite(initial$predicted))
    if (length(bad)) {
      stop("The initial orientation has no finite reference support for rows: ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    stop("The initial mean squared deviation is not finite. Rescale the measurements.",
         call. = FALSE)
  }
  best_offset <- c(0, 0, 0)
  best_metric <- initial$metric
  starts <- .orientation_starts(n_starts)
  attempts <- vector("list", n_starts)
  # A constant rescaling keeps optim's numerical tolerance independent of DN units.
  cost_scale <- if (initial$metric > 0) initial$metric else 1
  for (i in seq_len(n_starts)) {
    start <- starts[i, ]
    start_metric <- evaluate(start)$metric
    attempt_metric <- start_metric
    if (start_metric < best_metric) {
      best_metric <- start_metric
      best_offset <- start
    }
    objective <- function(offset) {
      if (any(!is.finite(offset)) || sqrt(sum(offset^2)) > 1) return(Inf)
      value <- evaluate(offset)$metric
      if (value < attempt_metric) attempt_metric <<- value
      if (value < best_metric) {
        best_metric <<- value
        best_offset <<- offset
      }
      value / cost_scale
    }
    convergence <- NA_integer_
    message <- "No finite reference support at this start."
    if (is.finite(start_metric)) {
      fit <- tryCatch(stats::optim(start, objective, method = method,
        control = list(maxit = max_iter, reltol = 1e-10,
                       parscale = rep(0.2, 3))), error = identity)
      if (inherits(fit, "error")) {
        message <- conditionMessage(fit)
      } else {
        convergence <- fit$convergence
        message <- if (is.null(fit$message)) "" else fit$message
      }
    }
    attempts[[i]] <- data.frame(east = start[1] * max_rotation,
      north = start[2] * max_rotation, up = start[3] * max_rotation,
      initial_metric = start_metric, metric = attempt_metric,
      convergence = convergence, message = message, row.names = NULL)
  }
  final <- evaluate(best_offset)
  orientation <- .orientation_from_frame(
    .orientation_rotation(best_offset * max_radians) %*% initial_frame)
  angles <- .angles_from_ray(final$external[, 1], final$external[, 2],
                             final$external[, 3], c("zenith", "azimuth"))
  samples <- sampling_directions
  samples$zenith <- angles[, 1]
  samples$azimuth <- angles[, 2]
  samples$reference_x <- final$x
  samples$reference_y_down <- final$y
  samples$predicted <- final$predicted
  samples$residual <- measured - final$predicted
  list(orientation = orientation, samples = samples, metric = final$metric,
       initial_metric = initial$metric,
       rotation_angle = sqrt(sum(best_offset^2)) * max_rotation,
       attempts = do.call(rbind, attempts))
}

.radiometric_orientation_frame <- function(orientation, name) {
  required <- c("sensor_rotation", "azimuth_rotation", "elevation")
  if (!is.list(orientation) || anyDuplicated(names(orientation)) ||
      !all(required %in% names(orientation))) {
    stop("`", name, "` must be a named list with sensor_rotation, ",
         "azimuth_rotation, and elevation.", call. = FALSE)
  }
  orientation <- orientation[required]
  do.call(.check_direction_orientation, orientation)
  basis <- data.frame(off_axis = c(pi/2, pi/2, 0),
                       internal_azimuth = c(pi/2, 0, 0))
  external <- do.call(external_directions_from_internal,
                      c(list(directions = basis), orientation))
  t(.orientation_unit_vectors(external, "zenith", "azimuth"))
}

.orientation_rotation <- function(vector) {
  angle <- sqrt(sum(vector^2))
  if (angle == 0) return(diag(3))
  axis <- vector / angle
  cross <- matrix(c(0, axis[3], -axis[2], -axis[3], 0, axis[1],
                    axis[2], -axis[1], 0), 3, 3)
  diag(3) + sin(angle) * cross + 2 * sin(angle/2)^2 * (cross %*% cross)
}

.orientation_starts <- function(n) {
  axes <- rbind(c(0, 0, 0), c(1, 0, 0), c(-1, 0, 0),
                c(0, 1, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1)) / 2
  if (n <= 7) return(axes[seq_len(n), , drop = FALSE])
  i <- seq_len(n - 7)
  z <- 1 - 2 * (i - 0.5) / length(i)
  phi <- i * pi * (3 - sqrt(5))
  rbind(axes, cbind(sqrt(1 - z^2) * cos(phi),
                    sqrt(1 - z^2) * sin(phi), z) / 2)
}

# Coordinates are measured from the upper-left image corner; centers are .5.
# Only neighbors with nonzero interpolation weight must have a finite value.
.orientation_bilinear <- function(pixels, x, y) {
  col <- x + 0.5
  row <- y + 0.5
  # Allow numerical round-trip error from radial inversion at pixel centers.
  tolerance <- sqrt(.Machine$double.eps) +
    32 * .Machine$double.eps * max(dim(pixels))
  near_col <- is.finite(col) & abs(col - round(col)) < tolerance
  near_row <- is.finite(row) & abs(row - round(row)) < tolerance
  col[near_col] <- round(col[near_col])
  row[near_row] <- round(row[near_row])
  valid <- is.finite(col) & is.finite(row) & col >= 1 & col <= ncol(pixels) &
    row >= 1 & row <= nrow(pixels)
  result <- rep(NA_real_, length(x))
  if (!any(valid)) return(result)
  col <- col[valid]
  row <- row[valid]
  c0 <- floor(col)
  r0 <- floor(row)
  c1 <- pmin(c0 + 1, ncol(pixels))
  r1 <- pmin(r0 + 1, nrow(pixels))
  dc <- col - c0
  dr <- row - r0
  weights <- cbind((1 - dc) * (1 - dr), dc * (1 - dr),
                    (1 - dc) * dr, dc * dr)
  values <- cbind(pixels[cbind(r0, c0)], pixels[cbind(r0, c1)],
                   pixels[cbind(r1, c0)], pixels[cbind(r1, c1)])
  baseline <- values[, 1]
  values[weights == 0] <- baseline[row(values)[weights == 0]]
  # Interpolate differences so that a constant reference stays exactly constant.
  upper <- baseline + dc * (values[, 2] - baseline)
  lower <- values[, 3] + dc * (values[, 4] - values[, 3])
  result[valid] <- upper + dr * (lower - upper)
  result
}
