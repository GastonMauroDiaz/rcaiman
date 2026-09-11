#' Fit a trend surface to sky digital numbers
#'
#' @description Fits a trend surface to sky digital numbers using
#' [spatial::surf.ls()] as the computational workhorse.
#'
#' @details This function models the variation in digital numbers across the sky
#' dome by fitting a polynomial surface in Cartesian space. It is intended to
#' capture smooth large-scale gradients and is more effective when called via
#' [apply_by_direction()] or the [panoramic()] pipeline.
#'
#' @param extrapolate logical vector of length one. If `TRUE`, predictions are
#'   extrapolated to the entire extent of `r`; otherwise, predictions are
#'   limited to the convex hull of the input sky points.
#' @param return character vector of length one. Either `"raster"` (default) to
#'   return the fitted surface as a raster, or `"all"` to also return the model
#'   and coefficient of determination.
#'
#' @inheritParams interpolate_planar
#' @inheritParams spatial::surf.ls
#'
#' @seealso [panoramic()] for an alternative pipeline involving a temporary
#'   shift of image plane
#'
#' @export
#'
#' @return
#' If `return = "raster"`, a numeric [terra::SpatRaster-class] with one layer
#' and the same geometry as `r`.
#'
#' If `return = "all"`, a list with named elements:
#' \describe{
#'   \item{`raster`}{[terra::SpatRaster-class] containing the fitted surface.}
#'   \item{`model`}{object of class `trls` returned by [spatial::surf.ls()].}
#'   \item{`r2`}{numeric value giving the coefficient of determination (R\eqn{^2}) of the fit.}
#' }
#'
#' @examples
#' r <- set_rcaiman_geometry(terra::rast(nrows = 5, ncols = 7))
#' sampling_coord <- coord_from_indices(
#'   col = c(1, 7, 1, 7, 2, 5, 3, 6),
#'   row = c(1, 1, 5, 5, 3, 2, 4, 4)
#' )
#' sampling_coord$dn <- c(2, 7, 4, 11, 3, 9, 6, 10)
#' fitted <- fit_trend_surface(
#'   sampling_coord, r = r, np = 2, extrapolate = TRUE, return = "all"
#' )
#' plot(fitted$raster)
#' fitted$r2
fit_trend_surface <- function(sampling_coord,
                              col_id = "dn",
                              r,
                              np = 6,
                              extrapolate = FALSE,
                              return = "raster") {

  if (!is.data.frame(sampling_coord)) {
    stop("`sampling_coord` must be a data frame.")
  }
  .assert_single_layer(r, allow_na = TRUE)
  .check_vector(np, "integerish", 1, sign = "positive")
  handling_col_id <- c(
    tryCatch(.check_vector(col_id, "numeric", sign = "any"),
             error = function(e) FALSE),
    tryCatch(.check_vector(col_id, "character"),
             error = function(e) FALSE)
  )
  if (!any(handling_col_id)) {
    stop("`col_id` must be the name or position of the column in `sampling_coord` containing the values to interpolate.")
  }
  if (is.numeric(col_id)) col_id <- names(sampling_coord)[col_id]
  required_cols <- c("x", "y_down", col_id)
  if (!all(required_cols %in% names(sampling_coord))) {
    stop(sprintf("`sampling_coord` must contain columns %s.",
                 paste(sprintf('"%s"', required_cols), collapse = ", ")))
  }
  .check_vector(extrapolate, "logical", 1)

  sampling_xy <- cbind(
    x = sampling_coord$x,
    y = y_from_y_down(sampling_coord$y_down, nrow(r))
  )
  fit <- spatial::surf.ls(x = sampling_xy[, "x"],
                          y = sampling_xy[, "y"],
                          z = sampling_coord[, col_id],
                          np)
  xl <- terra::xmin(r)
  xu <- terra::xmax(r)
  yl <- terra::ymin(r)
  yu <- terra::ymax(r)

  out <- spatial::trmat(fit, xl, xu, yl, yu, max(ncol(r), nrow(r)))
  out <- terra::rast(out$z) %>% t %>% flip
  terra::crs(out) <- terra::crs(r)
  terra::ext(out) <- terra::ext(r)
  out <- terra::resample(out, r)

  observed <- sampling_coord[, col_id]
  fitted <- terra::extract(out, sampling_xy)[,]
  ss_res <- sum((observed - fitted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  r2 <- 1 - ss_res / ss_tot

  if (!extrapolate) {
    hull_xy <- sampling_xy[grDevices::chull(sampling_xy),]
    v <- terra::vect(list(hull_xy), type = "polygon", crs = terra::crs(r))
    v <- terra::rasterize(v, r) %>% is.na()
    out[v] <- NA
  }
  if (return == "raster") {
    return(out)
  } else {
    return(list(raster = out, model = fit, r2 = r2))
  }

}
