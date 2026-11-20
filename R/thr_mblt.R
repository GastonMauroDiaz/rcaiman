#' Compute model-based thresholds
#'
#' @description
#' Compute threshold values from background digital numbers (DN) using
#' Equation 1 in \insertCite{Diaz2018;textual}{rcaiman}, a linear function
#' whose slope can be weighted.
#'
#' @details
#' For a Nikon Coolpix 5700 in JPEG mode the optimal threshold value is
#'  \eqn{\mathrm{OTV} = a + b \times W},
#' with \eqn{a = -7.8 / 255}, \eqn{b = 0.98},
#' and \eqn{W = 0.5}.
#'
#' The model was derived from canopy targets (perforated, rigid, dark
#' surfaces) backlit under homogeneous illumination and photographed with a
#' Nikon Coolpix 5700 in JPEG mode. Images were gamma-back-corrected with
#' a default gamma of 2.2 (see [invert_gamma_correction()]). Results showed that the optimal
#' threshold is linearly related to the background DN (see Figures 1 and 7
#' in \insertCite{Diaz2018;textual}{rcaiman}). This shifted the goal from
#' estimating an optimal threshold \insertCite{Song2014;textual}{rcaiman} to
#' estimating the background DN as if the canopy were absent, as proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. Currect versions of the method
#' developed with that intend can be found in [fit_trend_surface],
#' [fit_cone_shaped_model()], and [sky_via_local_models()].
#'
#' To apply the weighting parameter (W) from the equation, supply `slope` as
#' \eqn{slope \times w}.
#'
#' @note
#' Users can find all the details needed to calibrate their own cameras in
#' \insertCite{Diaz2018;textual}{rcaiman}, but they should acquire raw files
#' instead of JPEG, in such case this calibration is not required.
#'
#' @param dn numeric vector or [terra::SpatRaster-class]. Background digital
#'   number. Values must be gamma back corrected with [invert_gamma_correction].
#' @param intercept,slope numeric vectors of length one. Linear coefficients.
#'
#' @return An object of the same class and dimensions as `dn`.
#'
#' @seealso [normalize_min_max()], [invert_gamma_correction()]
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' thr_mblt(invert_gamma_correction(125))
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_jpeg.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' summary(my_entry)
#' jpeg_geometry <- get_geometry_spec(my_entry,  "jpeg", "rgb", "simple_method_rescaled")
#' jpeg_radiometry <- get_radiometry_spec(jpeg_geometry, "wagner2001_method")
#'
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path)
#' g <- set_internal_geometry(my_entry$jpeg$dim, jpeg_geometry)
#' ff <- flat_field_image(g$off_axis,
#'                        "simple_polynomial",
#'                        jpeg_radiometry$parameters$`0`$Blue)
#' plot(ff)
#' r <- apply_flat_field(caim$Blue, ff)
#' plot(r)
#' r <- invert_gamma_correction(r, 2.2)
#' plot(r)
#'
#' d <- direction_image(g)
#' seg <- beckers_segmentation(d, 20)
#' bin <- binarize_by_region(r, seg, "thr_two_corner")
#' plot(bin)
#' seg <- healpix_segmentation(d, 32)
#' seg[!bin] <- 0
#' w <- area_weight_image(g, jpeg_geometry)
#' features <- compute_features(distance_to_black(bin), seg, d, w)
#' sky_points <- features[features$effective_pixel_area > 3, ]
#' sky_points <- internal_directions_from_external(sky_points) %>%
#'   coord_from_internal_directions(., jpeg_geometry) %>% add_indices()
#' sky_points <- extract_dn(r, sky_points)
#' names(sky_points)
#' sky <- sky_via_local_models(d, sky_points, spacing_deg = 20, "Blue", conicity = 0.1)
#' plot(sky)
#'
#' thr <- thr_mblt(sky)
#' plot(binarize_with_thr(r, thr))
#' }
thr_mblt <- function (dn, intercept = -7.8 / 255, slope = 0.98 * 0.5) {
  handling_dn <- c(
    tryCatch(.check_vector(dn, "numeric", sign = "any"),
             error = function(e) FALSE),
    tryCatch(.assert_single_layer(dn),
             error = function(e) FALSE)
  )
  if (!any(handling_dn)) {
    stop("`dn` must be either numeric vector or SpatRaster class.")
  }
  .check_vector(intercept, "numeric", 1, sign = "any")
  .check_vector(slope, "numeric", 1, sign = "positive")

  thr <- intercept  + slope * dn
  thr[thr < 0] <- 0
  thr[is.na(thr)] <- 0

  thr
}
