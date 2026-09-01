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
#' @param return 	character vector of length one. Either "raster" (default) to
#'   return the fitted surface as a raster, or "all" to also return the model
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
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- extract_dn(r, sky_points)
#'
#' sky_s <- fit_trend_surface(sky_points, col_id = 3, r, np = 4,
#'                            extrapolate = TRUE)
#' plot(sky_s)
#' binarize_with_thr(r/sky_s, 0.5) %>% plot()
#'
#' sky_s <- fit_trend_surface(sky_points, col_id = 3, r, np = 6,
#'                            extrapolate = FALSE)
#' plot(sky_s)
#' binarize_with_thr(r/sky_s, 0.9) %>% plot()
#' }
fit_trend_surface <- function(sky_points,
                              col_id = "dn",
                              r,
                              np = 6,
                              extrapolate = FALSE,
                              return = "raster") {

  if (!is.data.frame(sky_points)) {
    stop("`sky_points` must be a data frame.")
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
    stop("`col_id` must be the name or position of the column in `sky_points` containing the values to interpolate.")
  }
  if (is.numeric(col_id)) col_id <- names(sky_points)[col_id]
  required_cols <- c("row", "col", col_id)
  if (!all(required_cols %in% names(sky_points))) {
    stop(sprintf("`sky_points` must contain columns %s.",
                 paste(sprintf('"%s"', required_cols), collapse = ", ")))
  }
  .check_vector(extrapolate, "logical", 1)

  xy <- terra::cellFromRowCol(r, sky_points$row, sky_points$col) %>%
    terra::xyFromCell(r, .)
  fit <- spatial::surf.ls(x = xy[, 1],
                          y = xy[, 2],
                          z = sky_points[, col_id],
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

  observed <- sky_points[, col_id]
  fitted <- terra::extract(out, xy)[,]
  ss_res <- sum((observed - fitted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  r2 <- 1 - ss_res / ss_tot

  if (!extrapolate) {
    xy <- xy[grDevices::chull(xy),]
    v <- terra::vect(list(xy), type = "polygon", crs = terra::crs(r))
    v <- terra::rasterize(v, r) %>% is.na()
    out[v] <- NA
  }
  if (return == "raster") {
    return(out)
  } else {
    return(list(raster = out, model = fit, r2 = r2))
  }

}
