#' Cone-shaped sky image
#'
#' @description
#' Evaluate parameters from [fit_cone_shaped_model()] at the external viewing
#' direction of each pixel.
#'
#' @details
#' The model formulas and their interpretation are described in
#' [fit_cone_shaped_model()]. The angles in `d` are converted from degrees to
#' radians before evaluation. Predictions use external directions regardless
#' of the image projection or acquisition orientation.
#'
#' Evaluation does not refit or normalize the model, and negative predictions
#' are retained. All available directions selected by `m` are evaluated,
#' including those outside the sampled angular range. Use the mask to limit
#' the region where predictions are appropriate.
#'
#' For `"zenith_only"`, the `azimuth` layer is ignored, including missing values.
#' For `"zenith_n_azimuth"`, both angles are required: a missing azimuth produces
#' `NA`, even at zenith. A finite azimuth at zenith is used as supplied; see
#' [fit_cone_shaped_model()] for the limitation of the sinusoidal terms there.
#'
#' @param d [terra::SpatRaster-class] from [direction_image()], with layers
#'   `zenith` and `azimuth` in degrees, in that order.
#' @param parameters named list returned by [fit_cone_shaped_model()], containing
#'   `method` and `coefficients`.
#' @param m optional logical [terra::SpatRaster-class] with one layer, no `NA`,
#'   and the same grid as `d`. Only cells where `m` is `TRUE` are evaluated.
#'   If `NULL`, all cells with the required directions are evaluated.
#'
#' @return numeric [terra::SpatRaster-class] with one layer named `prediction`,
#'   on the same grid as `d`. Cells outside the mask or missing a required angle
#'   are `NA`. Non-finite predictions at evaluated cells raise an error.
#'
#' @seealso [fit_cone_shaped_model()], [direction_image()]
#'
#' @note
#' This method was presented in \insertCite{Diaz2018;textual}{rcaiman}.
#'
#' @references \insertAllCited{}
#'
#' @export
#' @examples
#' # See fit_cone_shaped_model() for sampling an image and fitting the model.
#' samples <- sky_grid_centers(15)
#' samples$dn <- 1000 + 200 * samples$zenith^2
#' parameters <- fit_cone_shaped_model(samples, "dn", method = "zenith_only")
#' geometry <- list(model = "radial_projection", parameters = equidistant_coef,
#'   center_of_vision = c(40, 30), normalization_radius = 30,
#'   normalization_angle = 90)
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' m <- !is.na(d$zenith) & d$zenith < 90
#' sky <- cone_shaped_image(d, parameters, m = m)
#' plot(sky)
cone_shaped_image <- function(d, parameters, m = NULL) {
  .check_cone_shaped_parameters(parameters)
  z <- .sky_segmentation_layer(d, "zenith")
  zenith <- .degree2radian(terra::values(z, mat = FALSE))
  selected <- !is.na(zenith)
  if (parameters$method == "zenith_n_azimuth") {
    a <- .sky_segmentation_layer(d, "azimuth")
    azimuth <- .degree2radian(terra::values(a, mat = FALSE))
    selected <- selected & !is.na(azimuth)
  }
  if (!is.null(m)) {
    .assert_logical_mask(m)
    .assert_same_geom(m, z)
    selected <- selected & terra::values(m, mat = FALSE)
  }
  b <- parameters$coefficients
  theta <- zenith[selected]
  predicted <- b[["(Intercept)"]] + b[["zenith"]] * theta +
    b[["I(zenith^2)"]] * theta^2
  if (parameters$method == "zenith_n_azimuth") {
    phi <- azimuth[selected]
    predicted <- predicted + b[["sin(azimuth)"]] * sin(phi) +
      b[["cos(azimuth)"]] * cos(phi)
  }
  if (any(!is.finite(predicted))) {
    stop("Model predictions must be finite.", call. = FALSE)
  }
  response <- rep(NA_real_, terra::ncell(z))
  response[selected] <- predicted
  out <- terra::rast(z)
  terra::values(out) <- response
  names(out) <- "prediction"
  out
}

.check_cone_shaped_parameters <- function(parameters) {
  if (!is.list(parameters) ||
      !all(c("method", "coefficients") %in% names(parameters)) ||
      anyDuplicated(names(parameters))) {
    stop("`parameters` must be a named list with `method` and `coefficients`.",
         call. = FALSE)
  }
  .assert_choice(parameters$method, c("zenith_only", "zenith_n_azimuth"))
  expected <- c("(Intercept)", "zenith", "I(zenith^2)")
  if (parameters$method == "zenith_n_azimuth") {
    expected <- c(expected, "sin(azimuth)", "cos(azimuth)")
  }
  .check_vector(parameters$coefficients, "numeric", length(expected))
  if (!setequal(names(parameters$coefficients), expected)) {
    stop("`coefficients` must retain the names returned by `fit_cone_shaped_model()`.",
         call. = FALSE)
  }
  invisible(TRUE)
}
