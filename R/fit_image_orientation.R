#' Fit image orientation from corresponding landmarks
#'
#' @description
#' Estimate the orientation of an acquisition by matching its internal viewing
#' directions to the external directions of the same landmarks in a reference.
#'
#' @param internal_directions data frame with numeric columns `off_axis` and
#'   `internal_azimuth`, in radians, as returned by [add_directions()] with
#'   `g`. These are the directions in the image to orient.
#' @param reference_directions data frame with numeric columns `zenith` and
#'   `azimuth`, in radians, as returned by [add_directions()] with `d` for the
#'   reference acquisition. Rows must correspond to the same landmarks and be
#'   in the same order as `internal_directions`.
#'
#' @details
#' The acquisition procedure must place the entrance pupils at the same
#' location, so that corresponding landmarks differ only through orientation.
#' The function assumes this condition; it does not estimate translation or
#' correct parallax. Errors in geometric calibration and landmark selection
#' also contribute to the residuals.
#'
#' The reference directions establish the external frame. Their orientation
#' must already have been declared with [direction_image()]. Correspondences
#' alone cannot establish geographic north or the vertical.
#'
#' Directions are converted to unit vectors. A singular value decomposition
#' finds the rotation minimizing the sum of squared Euclidean distances
#' between paired vectors, with equal weight for every pair. Reflections are
#' excluded. This objective uses chord distances on the unit sphere; the
#' reported residuals are angular distances. No landmarks are removed and no
#' image resampling is involved.
#'
#' At least two pairs with nonparallel directions are required. Several
#' landmarks spread across the shared field of view allow the fit to be
#' assessed and avoid the instability of closely clustered directions. Remove
#' incomplete pairs from both tables together before fitting. Azimuth may be
#' `NA` at an angle of exactly 0 or pi, where the direction is defined without it.
#' Additional columns are ignored; pairs are matched by row, not by an ID.
#'
#' The fitted rotation is expressed using the convention of [direction_image()].
#' Both rotation arguments are in `[0, 360)` degrees, and elevation is in
#' `[-90, 90]` degrees. At a vertical optical axis, the two rotations are not
#' separately identifiable: the function sets `sensor_rotation = 0` and
#' expresses the sensor orientation through `azimuth_rotation`. Near vertical,
#' large changes in individual arguments can describe small changes in the
#' resulting orientation.
#'
#' The fitted orientation can serve as a starting point for visual refinement.
#' See "Describing an imperfectly levelled acquisition" in [direction_image()]
#' for a workflow using reprojection to a common destination and comparison
#' with flicker in rcaiman-projects.
#'
#' @return A plain list containing:
#'   \describe{
#'     \item{`orientation`}{a named list with `sensor_rotation`,
#'       `azimuth_rotation`, and `elevation`, in degrees, for [direction_image()].}
#'     \item{`residuals`}{a data frame in input row order, with fitted external
#'       `zenith` and `azimuth`, and `angular_error` relative to the reference.
#'       All three columns are in radians. Fitted azimuth is `NA` at zenith
#'       and nadir.}
#'     \item{`summary`}{a named numeric vector with `mean`, `rmse`, and `max`
#'       angular error, in radians.}
#'   }
#'
#' @seealso [add_directions()], [direction_image()], [reproject_image()],
#'   [optimize_orientation_by_radiometry()] for alignment against a sky
#'   reference using radiometric agreement instead of paired landmarks.
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection", parameters = equidistant_coef,
#'   center_of_vision = c(50, 50), normalization_radius = 50,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(100, 100), geometry)
#' # A known orientation supplies reference directions for this example.
#' d_reference <- direction_image(g, sensor_rotation = 90,
#'                                azimuth_rotation = 90, elevation = 87)
#' landmark_indices <- data.frame(col = c(25, 70, 60, 40),
#'                                row = c(35, 45, 75, 65))
#' internal <- add_directions(landmark_indices, g = g)
#' reference <- add_directions(landmark_indices, d = d_reference)
#' fit <- fit_image_orientation(internal, reference)
#' fit$orientation
#' fit$residuals$angular_error * 180 / pi
#' d <- do.call(direction_image, c(list(g = g), fit$orientation))
#'
#' # With two acquisitions, obtain each table from its own landmark indices
#' # and geometry. Keep corresponding landmarks in the same row order.
#' # internal <- add_directions(target_indices, g = target_g)
#' # reference <- add_directions(reference_indices, d = reference_d)
#' # fit <- fit_image_orientation(internal, reference)
#' # d <- do.call(direction_image,
#' #              c(list(g = target_g), fit$orientation))
fit_image_orientation <- function(internal_directions, reference_directions) {
  .check_columns(internal_directions,
                 required_columns = c("off_axis", "internal_azimuth"))
  .check_columns(reference_directions,
                 required_columns = c("zenith", "azimuth"))
  if (nrow(internal_directions) != nrow(reference_directions)) {
    stop("The direction tables must have the same number of rows.", call. = FALSE)
  }
  if (nrow(internal_directions) < 2L) {
    stop("At least two pairs of landmarks are required.", call. = FALSE)
  }
  internal <- .orientation_unit_vectors(internal_directions,
                                        "off_axis", "internal_azimuth")
  # Use the external frame of a default upward acquisition as the source
  # basis: sensor right points west, sensor top north, optical axis upward.
  internal[, 1] <- -internal[, 1]
  reference <- .orientation_unit_vectors(reference_directions,
                                         "zenith", "azimuth")
  decomposition <- svd(crossprod(internal, reference))
  handedness <- if (det(decomposition$v %*% t(decomposition$u)) < 0) -1 else 1
  tolerance <- 64 * .Machine$double.eps * decomposition$d[1]
  if (decomposition$d[2] <= tolerance ||
      (handedness < 0 &&
       decomposition$d[2] - decomposition$d[3] <= tolerance)) {
    stop("Landmark pairs do not determine a unique rotation. ",
         "Use nonparallel, corresponding directions.", call. = FALSE)
  }
  rotation <- decomposition$v %*% diag(c(1, 1, handedness)) %*%
    t(decomposition$u)
  # Columns now describe sensor right, sensor top, and optical axis in ENU.
  frame <- rotation %*% diag(c(-1, 1, 1))
  orientation <- .orientation_from_frame(frame)

  fitted <- internal %*% t(rotation)
  transverse <- sqrt(rowSums(fitted[, 1:2, drop = FALSE]^2))
  zenith <- atan2(transverse, fitted[, 3])
  azimuth <- atan2(fitted[, 1], fitted[, 2]) %% (2 * pi)
  vertical <- transverse <= 16 * .Machine$double.eps
  zenith[vertical] <- ifelse(fitted[vertical, 3] >= 0, 0, pi)
  azimuth[vertical] <- NA_real_
  azimuth[!is.na(azimuth) & azimuth >= 2 * pi] <- 0
  # atan2 of cross-product length and dot product is stable near zero error.
  cross <- cbind(fitted[, 2] * reference[, 3] - fitted[, 3] * reference[, 2],
                 fitted[, 3] * reference[, 1] - fitted[, 1] * reference[, 3],
                 fitted[, 1] * reference[, 2] - fitted[, 2] * reference[, 1])
  error <- atan2(sqrt(rowSums(cross^2)), rowSums(fitted * reference))
  list(orientation = orientation,
       residuals = data.frame(zenith = zenith, azimuth = azimuth,
                              angular_error = error),
       summary = c(mean = mean(error), rmse = sqrt(mean(error^2)),
                   max = max(error)))
}

# Shared conversion of sensor-right, sensor-top, optical-axis columns in ENU.
.orientation_from_frame <- function(frame) {
  axis <- frame[, 3]
  horizontal <- sqrt(sum(axis[1:2]^2))
  elevation <- atan2(axis[3], horizontal)
  if (horizontal <= 16 * .Machine$double.eps) {
    elevation <- sign(axis[3]) * pi / 2
    sensor_rotation <- 0
    azimuth_rotation <- atan2(sign(axis[3]) * frame[1, 2],
                              sign(axis[3]) * frame[2, 2])
  } else {
    azimuth_rotation <- atan2(-axis[1], -axis[2])
    right <- c(-cos(azimuth_rotation), sin(azimuth_rotation), 0)
    top <- c(sin(elevation) * sin(azimuth_rotation),
             sin(elevation) * cos(azimuth_rotation), cos(elevation))
    sensor_rotation <- atan2(sum(frame[, 2] * right), sum(frame[, 2] * top))
  }
  angles <- (c(sensor_rotation, azimuth_rotation) * 180 / pi) %% 360
  angles[angles >= 360] <- 0
  list(sensor_rotation = angles[1], azimuth_rotation = angles[2],
       elevation = elevation * 180 / pi)
}

.orientation_unit_vectors <- function(directions, polar_name, azimuth_name) {
  polar <- directions[[polar_name]]
  azimuth <- directions[[azimuth_name]]
  .check_vector(polar, "numeric", name = polar_name)
  .check_vector(azimuth, "numeric", allow_na = TRUE, name = azimuth_name)
  if (any(polar < 0 | polar > pi) ||
      any(azimuth < 0 | azimuth >= 2 * pi, na.rm = TRUE)) {
    stop("Directions must be in radians: polar angles in [0, pi] and ",
         "azimuths in [0, 2 * pi).", call. = FALSE)
  }
  axial <- polar %in% c(0, pi)
  if (any(is.na(azimuth) & !axial)) {
    stop("Azimuth may be NA only for a polar angle of 0 or pi. ",
         "Remove incomplete landmark pairs from both tables.", call. = FALSE)
  }
  azimuth[axial] <- 0
  vectors <- cbind(sin(polar) * sin(azimuth), sin(polar) * cos(azimuth), cos(polar))
  vectors[axial, 1:2] <- 0
  vectors
}
