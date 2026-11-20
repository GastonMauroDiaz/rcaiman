#' Compute local robust coefficient of variation image
#'
#' @description
#' Compute a raster of local robust coefficient of variation (CV) using a
#' fixed \eqn{3 \times 3} neighbourhood.
#'
#' @details
#' For each pixel, a local neighbourhood is extracted and a robust coefficient
#' of variation is computed as:
#'
#' \deqn{
#' \mathrm{CV}_{robust}
#' =
#' \frac{\mathrm{MAD}}{\mathrm{median}}
#' }
#'
#' where MAD is the median absolute deviation:
#'
#' \deqn{
#' \mathrm{MAD}
#' =
#' 1.4826
#' \times
#' \mathrm{median}
#' \left(
#' |x_i-\mathrm{median}(x)|
#' \right)
#' }
#'
#' Compared with mean-based coefficients of variation, this estimator is less
#' sensitive to outliers and mixed pixels. In hemispherical photography it may
#' help identify local radiometric instability caused by canopy edges,
#' sunlit elements, or abrupt local changes.
#'
#' Computation uses [terra::focalCpp()] with a C++ implementation optimized for
#' fixed \eqn{3 \times 3} windows.
#'
#' To obtain CV values at selected cells, compute the image and sample it
#' without applying another averaging window:
#'
#' ```
#' cv <- cv_image(r)
#' extract_dn(cv, sampling_indices, use_window = FALSE)
#' ```
#'
#' @param r Single-layer numeric [terra::SpatRaster-class].
#'
#' @return Single-layer numeric [terra::SpatRaster-class].
#'
#' @examples
#' r <- set_rcaiman_geometry(terra::rast(nrows = 10, ncols = 10, nlyrs = 1))
#' values(r) <- runif(ncell(r))
#' cv_image(r)
#'
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
#' ff <- flat_field_image(z, "simple_polynomial", rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' plot(r)
#'
#' cv <- cv_image()
#' plot(cv)
#' }
#' @export
cv_image <- function(r) {

  .assert_single_layer(r, allow_na = FALSE)

  cv <- terra::focalCpp(
    r,
    w = 3,
    fun = robust_cv
  )
  cv[is.na(cv)] <- 0
  names(cv) <- "cv"
  cv
}
