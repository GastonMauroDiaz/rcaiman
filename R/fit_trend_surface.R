
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
  out <- terra::rast(out$z)
  crs(out) <- crs(x)
  terra::ext(out) <- terra::ext(x)
  out <- terra::resample(out, x)

  list(image = out, fit = fit)
}



#' Fit a trend surface to sky digital numbers
#'
#' Fit a trend surface using spatial::surf.ls as workhorse function.
#'
#' This function is meant to be used after \code{\link{fit_coneshaped_model}}.
#'
#' A short explanation of this function can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}, after the explanation of the
#' \code{\link{fit_coneshaped_model}}.
#'
#' The argument \code{fact} is passed to \code{\link[terra]{aggregate}}. That
#' argument allows to control the scale at which the fitting is performed. In
#' general, a coarse scale lead to best generalization. The function used for
#' aggregation is \code{\link[stats]{quantile}}, to which the argument
#' \code{prob} is passed.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018}{rcaiman}.
#'
#' @inheritParams stats::quantile
#' @inheritParams ootb_mblt
#' @param bin \linkS4class{SpatRaster}. A working binarized image. This should
#'   be a preliminary binarization of \code{r}.
#' @param filling_source \linkS4class{SpatRaster}. Default is \code{NULL}.
#'   Above-canopy photograph. This image should contain pixels with sky DN
#'   values and \code{NA} in all the other pixels. A photograph taken
#'   immediately after or before taking \code{r} under the open sky with the
#'   same equipment and configuration is a very good option. The ideal option is
#'   one taken at the same time and place but above the canopy. The orientation
#'   relative to the North must be the same than for \code{r}.
#' @param prob Logical vector of length one. Probability for
#'   \code{\link[stats]{quantile}} calculation. See reference
#'   \insertCite{Diaz2018;textual}{rcaiman}.
#' @param m \linkS4class{SpatRaster}. A mask. Usually, the result of a call to
#'   \code{\link{mask_hs}}.
#' @inheritParams terra::aggregate
#' @inheritParams spatial::surf.ls
#'
#' @return A list with an object of class \linkS4class{SpatRaster} and of class
#'   \code{trls} (see \code{\link[spatial]{surf.ls}}).
#' @export
#'
#' @family MBLT functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' g <- sky_grid_segmentation(z, a, 10)
#' bin <- find_sky_pixels(blue, z, a)
#' sky_points <- extract_sky_points(blue, bin, g)
#' zenith_dn <- extract_zenith_dn(blue, z, a, sky_points)
#' rl_cs_fun <- fit_coneshaped_model(zenith_dn$sky_points)
#' sky_cs <- rl_cs_fun$rl_cs_fun(a, z) * zenith_dn$zenith_dn
#' m <- mask_hs(z, 0, 80)
#' sky <- fit_trend_surface(blue, bin, m, filling_source = sky_cs)
#' plot(sky$image)
#' }
fit_trend_surface <- function(r,
                              bin,
                              m = NULL,
                              filling_source = NULL,
                              prob = 0.95,
                              fact = 5,
                              np = 6) {
  .is_single_layer_raster(r, "r")
  .is_single_layer_raster(bin, "bin")
  .is_logic_and_NA_free(bin, "bin")
  terra::compareGeom(bin, r)
  if (!is.null(m)) {
    .is_single_layer_raster(m, "m")
    .is_logic_and_NA_free(m)
    terra::compareGeom(r, m)
    r[!m] <- NA
  }
  if (!is.null(filling_source)) {
    .is_single_layer_raster(filling_source, "filling_source")
    terra::compareGeom(r, filling_source)
  }
  stopifnot(length(prob) == 1)
  stopifnot(length(fact) == 1)
  stopifnot(length(np) == 1)

  fun <- function(x, ...) quantile(x, prob, na.rm = TRUE)

  blue <- r
  blue[!bin] <- NA
  if (fact > 1) blue <- terra::aggregate(blue, fact, fun, na.rm = TRUE)
  if (!is.null(filling_source)) {
    terra::compareGeom(bin, filling_source)
    if (fact > 1) {
      filling_source <- terra::aggregate(filling_source,
                                         fact, mean, na.rm = TRUE)
    }
    blue <- terra::cover(blue, filling_source)
  }

  surf <- .fit_trend_surface(blue, np = np)

  if (fact > 1) surf$image <- terra::resample(surf$image, r)
  if (!is.null(m)) surf$image[!m] <- NA

  surf
}
