
.fit_trend_surface <- function(x, np) {

  tmp <- data.frame(x = xFromCell(x, 1:ncell(x)),
                    y = yFromCell(x, 1:ncell(x)),
                    z = values(x)) %>%
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
#' This function is meant to be used before \code{\link{model_sky_dn}}.
#'
#' A short explanation of this function can be found on
#' \insertCite{Diaz2018;textual}{rcaiman}, under the heading \emph{Estimation of
#' the sky DN as a previous step for our method}, after the explanation of the
#' \code{\link{model_sky_dn}}.
#'
#' The argument \code{fact} is passed to the \code{\link[raster]{aggregate}}.
#' Whit this parameter, it is possible to control the scale at which the fitting
#' is performed. In general, a coarse scale lead to best generalization.
#'
#' If a stake of above canopy references is provided as filling source, the
#' function can select the one with less root mean square by comparing them on
#' the areas were sky DN are available in both \code{x} and
#' \code{filling_source}.
#'
#' @inheritParams model_sky_dn
#' @param mask raster
#' @inheritParams raster::aggregate
#' @inheritParams stats::quantile
#' @inheritParams spatial::surf.ls
#'
#' @return \code{\linkS4class{RasterLayer}}
#' @export
#'
#' @seealso mblt functions
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' a <- 10
fit_surface_to_sky_dn <- function(x,
                                  z,
                                  mask,
                                  bin,
                                  filling_source = NULL,
                                  prob = 0.95,
                                  fact = 5,
                                  np = 6) {


  compareRaster(bin, x)
  compareRaster(bin, mask)

  bin[!mask] <- NA

  fun <- function(x, ...) quantile(x, prob, na.rm = TRUE)

  Blue <- blue <- x * 255
  Blue[!bin] <- NA
  Blue[!mask] <- NA
  if (fact > 1) Blue <- raster::aggregate(Blue, fact, fun, na.rm = TRUE)

  if (!is.null(filling_source)){
    compareRaster(bin, filling_source)

    filling_source[!mask] <- NA
    filling_source <- filling_source * 255
    if (fact > 1) {
      filling_source <- raster::aggregate(filling_source,
                                          fact,
                                          mean,
                                          na.rm = TRUE)
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
  model <- r$fit
  r <- .filter_values(r$image)

  if (fact > 1) r <- resample(r, blue)

  # # 10/9/2019
  # # A new method to increase the robustness of the fit.
  # # It uses a plane to model the DNs near the horizon.
  #
  # max_angle <- round(max(z[mask]))
  # aux_mask <-  mask_image(z, zlim = round(c(max_angle - 1,
  #                                           max_angle)))
  #
  # aux_r <- r
  # aux_r[!aux_mask] <- NA
  # aux_r <- median(aux_r[], na.rm = TRUE) - IQR(aux_r[], na.rm = TRUE) %>%
  #   .filter_values(aux_r, ., NULL)
  #
  #
  # plane <- .fit_trend_surface(aux_r, 1)$image
  # plane[plane < 1] <- 1
  # plane[plane > 255] <- 255
  #
  # if (max_angle < 80) {
  #   aux_mask <- mask_image(z, zlim = c(0, max_angle + 10))
  # } else {
  #   aux_mask <- mask_image(z)
  # }
  # plane[aux_mask] <- NA
  #
  # if (fact > 1) plane <- raster::aggregate(plane, fact, fun, na.rm = TRUE)
  #
  # Blue <- cover(Blue, plane)
  #
  # r <- .fit_trend_surface(Blue, np = np)$image
  # r <- .filter_values(r)
  #
  # if (fact > 1) r <- resample(r, blue)

  list(image = r / 255, model = model)
}
