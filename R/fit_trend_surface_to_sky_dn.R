
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

.filter_values <- function(r, mn = 0, mx = 255) {
  if (!is.null(mn)) r[r < mn] <- NA
  if (!is.null(mx)) r[r > mx] <- NA
  r
}

#' Fit a trend surface to sky digital numbers
#'
#' Fit a trend surface using spatial::surf.ls as workhorse function.
#'
#' This function is meant to be used after \code{\link{model_sky_dn}}.
#'
#' A short explanation of this function can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}, after the explanation of the
#' \code{\link{model_sky_dn}}.
#'
#' The argument \code{fact} is passed to \code{\link[raster]{aggregate}}. That
#' argument allows to control the scale at which the fitting is performed. In
#' general, a coarse scale lead to best generalization.
#'
#' If a raster stake of above canopy references is provided as filling source,
#' \code{fit_trend_surface_to_sky_dn()} function can select the one with less root
#' mean square error by analyzing the areas were sky DN are available in both
#' \code{x} and \code{filling_source}. If no above canopy reference is available
#' the result from a call to \code{\link{model_sky_dn}} can be provided as
#' filling source. However, the function can work without a filling source,
#' which is the default setup.
#'
#' @inheritParams model_sky_dn
#' @param mask raster
#' @inheritParams raster::aggregate
#' @inheritParams stats::quantile
#' @inheritParams spatial::surf.ls
#'
#' @return A list with an object of class \linkS4class{RasterLayer} and of
#'   class "trls" (see \code{\link[spatial]{surf.ls}}).
#' @export
#'
#' @family mblt functions
#'
#' @references
#' \insertRef{Diaz2018}{rcaiman}
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
#' sky <- model_sky_dn(blue, z, a, bin, parallel = FALSE)
#' aux_m <- mask_image(z, zlim = c(0, 20))
#' sky$image[aux_m] <- NA
#' m <- mask_image(z, zlim = c(0, 70))
#' sky <- fit_trend_surface_to_sky_dn(blue, z, m, bin,
#'   filling_source = sky$image
#' )
#' }
fit_trend_surface_to_sky_dn <- function(r,
                                        z,
                                        mask,
                                        bin,
                                        filling_source = NULL,
                                        prob = 0.95,
                                        fact = 5,
                                        np = 6) {
  .check_if_r_was_normalized(r)

  compareRaster(bin, r)
  compareRaster(bin, mask)

  bin[!mask] <- NA

  fun <- function(x, ...) quantile(x, prob, na.rm = TRUE)

  Blue <- blue <- r * 255
  rm(r)
  Blue[!bin] <- NA
  Blue[!mask] <- NA
  if (fact > 1) Blue <- raster::aggregate(Blue, fact, fun, na.rm = TRUE)

  if (!is.null(filling_source)) {
    compareRaster(bin, filling_source)

    filling_source[!mask] <- NA
    filling_source <- filling_source * 255
    if (fact > 1) {
      filling_source <- raster::aggregate(filling_source,
        fact,
        mean,
        na.rm = TRUE
      )
    }
    Filling <- Map(
      function(x) raster::subset(filling_source, x),
      1:nlayers(filling_source)
    )

    # correct bias
    .findBias <- function(x, y) {
      m <- mask_image(z, zlim = c(30, 60))
      mean(x[m], na.rm = TRUE) - mean(y[m], na.rm = TRUE)
    }
    Filling <- Map(function(x) x - .findBias(x, Blue), Filling)

    # select the best filling data
    rmse <- function(error) sqrt(mean(error^2, na.rm = TRUE))
    RMSE <- Map(function(x) rmse((x - Blue)[]), Filling)
    index <- which.min(unlist(RMSE))
    Filling <- Filling[[index]]
    RMSE <- RMSE[[index]]

    # force corrected values
    Filling[Filling > 255] <- 255
    Filling[Filling < 0] <- 0

    # fill
    foo <- raster::sampleRegular(Filling, ncell(Filling) * 0.7,
      cells = TRUE
    )
    Filling[foo[, 1]] <- NA

    Blue <- cover(Blue, Filling)
  }

  r <- .fit_trend_surface(Blue, np = np)
  fit <- r$fit
  r <- .filter_values(r$image)

  if (fact > 1) r <- resample(r, blue)

  r[is.na(z)] <- NA

  list(image = r / 255, fit = fit)
}
