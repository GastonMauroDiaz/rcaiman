
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
#' @inheritParams stats::quantile
#' @inheritParams fit_cone_shaped_model
#' @param m \linkS4class{RasterLayer}. A mask. Usually, the result of a call to
#'   \code{\link{mask_image}}.
#' @inheritParams raster::aggregate
#' @inheritParams spatial::surf.ls
#'
#' @return A list with an object of class \linkS4class{RasterLayer} and of class
#'   "trls" (see \code{\link[spatial]{surf.ls}}).
#' @export
#'
#' @family mblt functions
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' \dontrun{
#' path <- getwd()
#' my_file <- paste0(path, "/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'   method = "auto", mode = "wb"
#' )
#' r <- read_caim(
#'   my_file,
#'   c(1280, 960) - 745,
#'   745 * 2,
#'   745 * 2
#' )
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' thr <- autothresholdr::auto_thresh(r$Blue[], "IsoData")
#' bin <- apply_thr(r$Blue, thr[1] * 1.25)
#' blue <- gbc(r$Blue)
#' sky <- fit_cone_shaped_model(blue, z, a, bin, parallel = FALSE)
#' m <- mask_image(z, zlim = c(0, 80))
#' sky <- fit_trend_surface(blue, m, bin,
#'   filling_source = sky$image
#' )
#' plot(sky$image)
#' }
fit_trend_surface <- function(r,
                              m,
                              bin,
                              filling_source = NULL,
                              prob = 0.95,
                              fact = 5,
                              np = 6) {
  stopifnot(class(r) == "RasterLayer")
  .check_if_r_was_normalized(r)
  compareRaster(bin, r)
  compareRaster(bin, m)

  r[!m] <- NA

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
  fit <- surf$fit

  surf <- surf$image
  surf[surf < 0] <- NA
  surf[surf > 1] <- NA

  if (fact > 1) surf <- resample(surf, r)
  surf[!m] <- NA

  list(image = surf , fit = fit)
}
