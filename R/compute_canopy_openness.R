#' Calculate canopy openness over a selected region
#'
#' @description
#' Calculate the solid-angle-weighted proportion of canopy gaps within the
#' selected observed coverage of the upper hemisphere.
#'
#' @param bin logical [terra::SpatRaster-class] with one layer. `TRUE` denotes
#'   a gap and `FALSE` vegetation. Missing values are excluded from the
#'   calculation. See [binarize_with_thr()].
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees and on the same grid as `bin`. Only
#'   `zenith` is used: angles from 0 to 90 deg select the upper hemisphere.
#'   Missing azimuths do not affect the calculation.
#' @param w numeric [terra::SpatRaster-class] with one layer and on the same grid
#'   as `bin`. Nonnegative, finite solid-angle weights, such as those returned
#'   by [area_weight_image()]. Relative weights are sufficient. Missing weights
#'   are excluded from the calculation.
#' @param m optional logical [terra::SpatRaster-class] with one layer, on the
#'   same grid as `bin`, and without missing values. `TRUE` selects pixels within
#'   the upper hemisphere. When `NULL`, all observed upper-hemisphere pixels
#'   with available `bin` and `w` values are used.
#'
#' @details
#' The solid-angle weighting approach follows the formulation implemented in
#' CIMES by Jean-Michel Walter \insertCite{Gonsamo2011}{rcaiman}. Here, weights
#' are assigned to individual pixels rather than angular segments and
#' normalized over the selected region, allowing canopy openness to be
#' estimated within partial hemispherical coverage.
#'
#' For the selected pixels \eqn{M}, openness is
#' \deqn{CO_M = \frac{\sum_{i \in M} b_i w_i}{\sum_{i \in M} w_i},}
#' where \eqn{b_i} is 1 for a gap and 0 for vegetation, and \eqn{w_i} is
#' proportional to the solid angle represented by that pixel. Multiplying all
#' weights by the same positive constant leaves the result unchanged. With
#' equal weights, the result is the ordinary gap fraction of the selected pixels.
#'
#' The denominator always refers to the selected valid coverage, whether or
#' not `m` is supplied. Missing observations are excluded from both sums;
#' the function does not estimate openness in unobserved directions. It stops
#' if the selection contains no positive weight.
#'
#' External directions in `d` delimit the upper hemisphere, independently of
#' acquisition orientation. Directions below the horizon or with missing
#' zenith angles are excluded, even when `m` selects them. In particular,
#' `!is.na(d$zenith)` alone does not identify the upper hemisphere.
#' Use `m` to further restrict the region, for example to zenith angles
#' from 0 to 70 deg.
#'
#' Weights describe internal geometry and do not change with acquisition
#' orientation. If inclining the system changes the observed portion of
#' the canopy, its openness can differ; unobserved directions are not
#' reconstructed or included in the denominator.
#'
#' No segmentation or reprojection is performed. With weights from
#' [area_weight_image()], the calculation uses the local solid-angle
#' approximation documented there and currently supports radial models.
#' Other valid solid-angle weights can also be supplied.
#'
#' @return Numeric vector of length one in \eqn{[0,1]}: the proportion of
#'   selected valid angular coverage occupied by gaps.
#' @seealso [area_weight_image()], [direction_image()]
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
#' d <- direction_image(g, elevation = 32.5)
#' w <- area_weight_image(g, geometry)
#'
#' # Illustrative binary pattern; in practice, classify the observed image.
#' bin <- g$internal_azimuth < 180
#' # The horizon is applied automatically, including in inclined acquisitions.
#' compute_canopy_openness(bin, d, w)
#'
#' # Restrict the calculation to a smaller external angular region.
#' m <- !is.na(d$zenith) & d$zenith <= 70
#' compute_canopy_openness(bin, d, w, m)
compute_canopy_openness <- function(bin, d, w, m = NULL) {
  .assert_single_layer(bin)
  z <- .sky_segmentation_layer(d, "zenith")
  .assert_same_geom(bin, z, y_name = "d")
  .assert_area_weight_image(w)
  .assert_same_geom(bin, w)
  gaps <- terra::values(bin, mat = FALSE)
  weights <- terra::values(w, mat = FALSE)
  .check_vector(gaps, "logical", allow_na = TRUE, name = "bin")
  zenith <- terra::values(z, mat = FALSE)
  selected <- !is.na(gaps) & !is.na(weights) &
    !is.na(zenith) & zenith <= 90
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(bin, m)
    selected <- selected & terra::values(m, mat = FALSE)
  }
  weights <- weights[selected]
  if (!length(weights) || !any(weights > 0)) {
    stop("The selected valid coverage must contain positive weights.",
         call. = FALSE)
  }
  # A common scale cancels in the ratio; rescale to avoid overflowing sums.
  weights <- weights / max(weights)
  sum(weights[gaps[selected]]) / sum(weights)
}
