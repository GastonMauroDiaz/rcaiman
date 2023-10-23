#' Fit a trend surface to sky digital numbers
#'
#' Fit a trend surface using [spatial::surf.ls()] as workhorse
#' function.
#'
#' This function is meant to be used after [fit_coneshaped_model()].
#'
#' This method was presented in \insertCite{Diaz2018;textual}{rcaiman}, under the
#' heading *Estimation of the sky DN as a previous step for our method*. If you
#' use this function in your research, please cite that paper in addition to
#' this package (`citation("rcaiman"`).
#'
#'
#' @inheritParams ootb_mblt
#' @param filling_source [SpatRaster-class]. An actual or reconstructed
#'   above-canopy image to complement the sky pixels detected through the gaps
#'   of `r`. If an incomplete above-canopy image is available, non-sky
#'   pixels should be turned `NA` or they will be considered as sky pixels
#'   erroneously. A photograph taken immediately after or before taking `r`
#'   under the open sky with the same equipment and configuration is a very good
#'   option but not recommended under fleeting clouds. The orientation relative
#'   to the North must be the same as for `r`. If it is set to `NULL`
#'   (default), only sky pixels from `r` will be used as input.
#' @inheritParams spatial::surf.ls
#'
#' @return A list with an object of class [SpatRaster-class] and of class
#'   `trls` (see [spatial::surf.ls()]).
#' @export
#'
#' @family Sky Reconstruction Functions
#' @seealso [thr_mblt()]
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' caim <- normalize(caim, 0, 20847, TRUE)
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' bin <- ootb_obia(caim, z, a, m, HSV(239, 0.85, 0.5), gamma = NULL)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_plant = 5)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' rl <- extract_rl(r, z, a, sky_points)
#'
#' model <- fit_coneshaped_model(rl$sky_points)
#' summary(model$model)
#' sky_cs <- model$fun(z, a)
#' persp(terra::aggregate(sky_cs, 10), theta = 90, phi = 45)
#'
#' sky_s <- fit_trend_surface(r, z, a, bin, sky_cs)
#' persp(terra::aggregate(sky_s$image, 10), theta = 90, phi = 45)
#'
#' #a quick pipeline
#' bin <- find_sky_pixels(r, z, a)
#' sky <- fit_trend_surface(r, z, a, bin)$image
#' sky <- fix_reconstructed_sky(sky, z, r, bin)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels_nonnull(r, sky, g)
#' plot(bin)
#' }
fit_trend_surface <- function(r,
                              z,
                              a,
                              bin,
                              filling_source = NULL,
                              np = 6) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  terra::compareGeom(bin, r)
  if (!is.null(filling_source)) {
    .is_single_layer_raster(filling_source, "filling_source")
    terra::compareGeom(r, filling_source)
  }
  stopifnot(length(np) == 1)

  r[!bin] <- NA
  if (!is.null(filling_source)) {
    terra::compareGeom(bin, filling_source)
    r <- terra::cover(r, filling_source)
  }

  .fit_trend_surface <- function(x, np) {
    tmp <- data.frame(
      terra::xFromCell(x, 1:ncell(x)),
      terra::yFromCell(x, 1:ncell(x)),
      terra::values(x)
    ) %>%
      .[!is.na(.[,3]), ]
    colnames(tmp) <- c("x", "y", names(x))
    fit <- spatial::surf.ls(x = tmp[, 1], y = tmp[, 2], z = tmp[, 3], np)
    xl <- xmin(x)
    xu <- xmax(x)
    yl <- ymin(x)
    yu <- ymax(x)

    out <- spatial::trmat(fit, xl, xu, yl, yu, ncol(x))
    out <- terra::rast(out$z) %>% t %>% flip
    crs(out) <- crs(x)
    terra::ext(out) <- terra::ext(x)
    out <- terra::resample(out, x)

    list(image = out, model = fit)
  }
  out <- .fit_trend_surface(r, np)
  out$image[is.na(z)] <- NA
  out
}

# pano = FALSE
# if (pano) {
#   r[!bin] <- filling_source[!bin]
#   angle_width <- 5
#   r <- fisheye_to_pano(r, z, a, angle_width = angle_width)
#   out <- .fit_trend_surface(r, np)
#   image <- out$image
#   .a <- .z <- r
#   terra::values(.a) <- .col(dim(r)[1:2])
#   terra::values(.z) <- .row(dim(r)[1:2])
#   g <- sky_grid_segmentation(z, a, angle_width)
#   cells <- 1:terra::ncell(image)
#   fun <- function(s, r) s * 1000 + r
#   .labels <- fun(.a[cells], .z[cells])
#   image <- terra::subst(g, .labels, image[cells])
#   out$image <- image
# } else {
#   out <- .fit_trend_surface(r, np)
# }
