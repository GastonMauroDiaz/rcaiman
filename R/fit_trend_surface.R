#' Fit a trend surface to sky digital numbers
#'
#' Fit a trend surface using \code{\link[spatial]{surf.ls}} as workhorse
#' function.
#'
#' This function is meant to be used after \code{\link{fit_coneshaped_model}}.
#'
#' A short explanation of this function can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}, after the explanation of the
#' \code{\link{fit_coneshaped_model}}.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018;textual}{rcaiman} in addition to this package.
#'
#' @inheritParams ootb_mblt
#' @param filling_source \linkS4class{SpatRaster}. An actual or reconstructed
#'   above-canopy image to complement the sky pixels detected through the gaps
#'   of \code{r}. If an incomplete above-canopy image is available, non-sky
#'   pixels should be turned \code{NA} or they will be considered as sky pixels
#'   erroneously. A photograph taken immediately after or before taking \code{r}
#'   under the open sky with the same equipment and configuration is a very good
#'   option but not recommended under fleeting clouds. The orientation relative
#'   to the North must be the same as for \code{r}. If it is set to \code{NULL}
#'   (default), only sky pixels from \code{r} will be used as input.
#' @inheritParams spatial::surf.ls
#'
#' @return A list with an object of class \linkS4class{SpatRaster} and of class
#'   \code{trls} (see \code{\link[spatial]{surf.ls}}).
#' @export
#'
#' @family Sky Reconstruction Functions
#' @seealso \code{\link{thr_image}}
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \donttest{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' r <- gbc(caim$Blue)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels(r, z, a)
#' sky_points <- extract_sky_points(r, bin, g)
#' sky_points <- extract_rl(r, z, a, sky_points, NULL)
#' model <- fit_coneshaped_model(sky_points$sky_points)
#' sky_cs <- model$fun(z, a)
#' m <- mask_hs(z, 0, 80)
#' sky <- fit_trend_surface(r, z, a, bin, filling_source = sky_cs)
#' plot(sky$image)
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

    list(image = out, fit = fit)
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
