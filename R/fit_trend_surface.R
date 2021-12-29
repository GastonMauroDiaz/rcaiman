
.fit_trend_surface <- function(x, np) {
  tmp <- data.frame(
    x = xFromCell(x, 1:ncell(x)),
    y = yFromCell(x, 1:ncell(x)),
    z = values(x)
  ) %>%
    .[!is.na(.$z), ]
  colnames(tmp) <- c("x", "y", names(x))

  fit <- spatial::surf.ls(x = tmp[, 1], y = tmp[, 2], z = tmp[, 3], np)
  xl <- xmin(x)
  xu <- xmax(x)
  yl <- ymin(x)
  yu <- ymax(x)

  out <- spatial::trmat(fit, xl, xu, yl, yu, ncol(x))
  out <- raster(out)
  out <- resample(out, x)

  list(image = out, fit = fit)
}



#' Fit a trend surface to sky digital numbers
#'
#' Fit a trend surface using spatial::surf.ls as workhorse function.
#'
#' This function is meant to be used after \code{\link{fit_cone_shaped_model}}.
#'
#' A short explanation of this function can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}, after the explanation of the
#' \code{\link{fit_cone_shaped_model}}.
#'
#' The argument \code{fact} is passed to \code{\link[raster]{aggregate}}. That
#' argument allows to control the scale at which the fitting is performed. In
#' general, a coarse scale lead to best generalization. The function used for
#' aggregation is \code{\link[stats]{quantile}}, to which the argument
#' \code{prob} is passed. Essentially, the aggregation step works as the one
#' from \code{\link{fit_cone_shaped_model}}, but it is made on the raster space
#' rather than on the hemispherical space.
#'
#' If you use this function in your research, please cite
#' \insertCite{Diaz2018}{rcaiman}.
#'
#' @inheritParams stats::quantile
#' @inheritParams fit_cone_shaped_model
#' @param m \linkS4class{RasterLayer}. A mask. Usually, the result of a call to
#'   \code{\link{mask_hs}}.
#' @inheritParams raster::aggregate
#' @inheritParams spatial::surf.ls
#'
#' @return A list with an object of class \linkS4class{RasterLayer} and of class
#'   \code{trls} (see \code{\link[spatial]{surf.ls}}).
#' @export
#'
#' @family MBLT functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/4_D_2_DSCN4502.JPG", package = "rcaiman")
#' r <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' bin <- find_sky_dns(blue, z, a)
#' sky <- fit_cone_shaped_model(blue, z, a, bin, parallel = FALSE)
#' m <- mask_hs(z, 0, 80)
#' sky <- fit_trend_surface(blue, bin, m, filling_source = sky$image)
#' plot(sky$image)
#' }
fit_trend_surface <- function(r,
                              bin,
                              m = NULL,
                              filling_source = NULL,
                              prob = 0.95,
                              fact = 5,
                              np = 6) {
  stopifnot(class(r) == "RasterLayer")
  compareRaster(bin, r)
  if (!is.null(m)) compareRaster(bin, m)

  if (!is.null(m)) r[!m] <- NA

  fun <- function(x, ...) quantile(x, prob, na.rm = TRUE)

  blue <- r
  blue[!bin] <- NA
  if (fact > 1) blue <- aggregate(blue, fact, fun, na.rm = TRUE)
  if (!is.null(filling_source)) {
    compareRaster(bin, filling_source)
    if (fact > 1) {
      filling_source <- aggregate(filling_source, fact, mean, na.rm = TRUE)
    }
    blue <- cover(blue, filling_source)
  }

  surf <- .fit_trend_surface(blue, np = np)

  if (fact > 1) surf$image <- resample(surf$image, r)
  if (!is.null(m)) surf$image[!m] <- NA

  surf
}
