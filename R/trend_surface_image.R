#' Trend surface image
#'
#' @description
#' Evaluate a fitted polynomial surface directly at image pixel centers.
#'
#' @details
#' Supply parameters from [fit_trend_surface()]. Predictions use the same
#' `x`, `y_down` coordinate frame as the observations, with the first pixel
#' center at `(0.5, 0.5)`. Evaluation does not refit the model, renormalize its
#' values, or resample an intermediate raster. Predictions may be negative
#' or zero.
#'
#' All unmasked cells are evaluated, including positions outside the sample
#' ranges. Use [convex_hull_mask()] or another mask to limit the result when
#' appropriate. Changing the grid dimensions does not rescale the model.
#' [crop_caim()] changes the coordinate origin: to retain correspondence,
#' predict in the original frame and crop the result together with the image.
#'
#' @section Radiometric preprocessing:
#' A fitted surface can represent the spatial sensitivity measured in a
#' flat-field calibration. Select observations for one band and aperture,
#' normalized to the response at the center of vision, on the optical axis.
#' Use `relative_response` for this instrumental response; it is distinct
#' from relative scene radiance.
#'
#' Fit that column with [fit_trend_surface()] and retain its parameter list.
#' Normalizing observations does not constrain the fitted surface to equal
#' one at the center of vision. Predictions are not renormalized.
#'
#' Evaluate the model on images corresponding to that calibration, using the
#' same sensor coordinate frame. Then apply the flat-field with
#' [apply_flat_field()], which checks positive, finite correction values.
#' Correct radiometry before reprojection. In a registry radiometry
#' specification, `model = "trend_surface_image"` identifies this general
#' evaluator; the specification supplies the radiometric meaning.
#'
#' The following example illustrates this application with synthetic data.
#'
#' @param r [terra::SpatRaster-class] with one layer and the unit-resolution
#'   grid established by [set_rcaiman_geometry()]. Its grid defines the output.
#'   If it contains values, `NA` cells mask the output; other values are ignored.
#'   A template without values is also accepted.
#' @param parameters named list returned by [fit_trend_surface()], containing
#'   `np`, `beta`, `rx`, and `ry`.
#' @param m optional logical [terra::SpatRaster-class] with one layer, no `NA`,
#'   and the same grid as `r`. Predictions are evaluated where `m` is `TRUE`
#'   and `r` is not `NA`. Other cells are `NA`. If `NULL`, only the template
#'   mask is used. See [convex_hull_mask()] for one way to construct `m`.
#'
#' @return numeric [terra::SpatRaster-class] with one layer named `prediction`
#'   and the same grid as `r`. Non-finite predictions in unmasked cells raise
#'   an error.
#'
#' @seealso [fit_trend_surface()], [convex_hull_mask()]
#' @export
#' @examples
#' # A normalized instrumental response is one application of the surface.
#' samples <- expand.grid(x = seq(0.5, 19.5, length.out = 8),
#'                        y_down = seq(0.5, 11.5, length.out = 8))
#' samples$relative_response <- 1 - 0.001 * (samples$x - 10)^2 -
#'   0.002 * (samples$y_down - 6)^2
#' parameters <- fit_trend_surface(samples, col_id = "relative_response", np = 2)
#' r <- set_rcaiman_geometry(terra::rast(nrows = 12, ncols = 20))
#' flat_field <- trend_surface_image(r, parameters)
#' # Simulate a uniform scene measured through this response.
#' r <- 100 * flat_field
#' corrected <- apply_flat_field(r, flat_field)
#' range(terra::values(corrected))
trend_surface_image <- function(r, parameters, m = NULL) {
  .check_image_template(r)
  .check_trend_surface_parameters(parameters)
  cells <- seq_len(terra::ncell(r))
  if (terra::hasValues(r)) {
    cells <- cells[!is.na(terra::values(r, mat = FALSE))]
  }
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, r)
    cells <- cells[terra::values(m, mat = FALSE)[cells]]
  }
  response <- rep(NA_real_, terra::ncell(r))
  if (length(cells)) {
    xy <- terra::xyFromCell(r, cells)
    fitted <- .predict_trend_surface(
      parameters, xy[, 1], terra::nrow(r) - xy[, 2]
    )
    if (any(!is.finite(fitted))) {
      stop("Surface predictions must be finite.",
           call. = FALSE)
    }
    response[cells] <- fitted
  }
  out <- terra::rast(r)
  terra::values(out) <- response
  names(out) <- "prediction"
  out
}

.check_image_template <- function(r) {
  if (!inherits(r, "SpatRaster")) {
    stop("`r` must be a SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(r) != 1) {
    stop("`r` must have a single layer.", call. = FALSE)
  }
  expected_extent <- c(0, terra::ncol(r), 0, terra::nrow(r))
  if (!isTRUE(all.equal(unname(as.vector(terra::ext(r))), expected_extent))) {
    stop("`r` must use the image grid from `set_rcaiman_geometry()`.",
         call. = FALSE)
  }
  invisible(TRUE)
}
