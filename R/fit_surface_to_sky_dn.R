
.fitTrendSurface <- function(x, sampleProportion, np) {
  for (i in seq(length = nlayers(x))) {
    if (nlayers(x) > 1) {
      tmp <- sampleRandom(raster::subset(x, i), ncell(x) * sampleProportion,
        sp = TRUE
      )
    } else {
      tmp <- sampleRandom(x, ncell(x) * sampleProportion, sp = TRUE)
    }

    tmp <- cbind(tmp@coords, tmp@data)

    fit <- spatial::surf.ls(x = tmp[, 1], y = tmp[, 2], z = tmp[, 3], np)
    xl <- xmin(x)
    xu <- xmax(x)
    yl <- ymin(x)
    yu <- ymax(x)

    out <- spatial::trmat(fit, xl, xu, yl, yu, ncol(x))
    out <- raster(out)
    out <- resample(out, x)
    if (nlayers(x) > 1) {
      x[[i]] <- out
    } else {
      x <- out
    }
  }
  list(x, fit)
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
fit_surface_to_sky_dn <- function(x,
                                  bin,
                                  mask,
                                  filling,
                                  prob = 0.95,
                                  fact = 5,
                                  np = 6) {


  compareRaster(bin, filling)
  compareRaster(bin, x)
  compareRaster(bin, mask)

  bin[!mask] <- NA
  filling[!mask] <- NA

  fun <- function(x, ...) quantile(x, prob, na.rm = TRUE)

  blue <- x

  Blue <- blue
  Blue[!bin] <- NA
  Blue[!mask] <- NA

  if (fact > 1) {
    Blue <- raster::aggregate(Blue, fact, fun, na.rm = TRUE)
    filling <- raster::aggregate(filling, fact, mean, na.rm = TRUE)
  }
  Filling <- Map(
    function(x) raster::subset(filling, x),
    1:nlayers(filling)
  )


  # correct bias in the data for filling
  .findBias <- function(x, y) {
    m <- mask_image(z, zlim = c(30, 60))
    mean(x[m], na.rm = TRUE) - mean(y[m], na.rm = TRUE)
  }
  Filling <- Map(function(x) x - .findBias(x, Blue), Filling)

  # select the better data for filling
  rmse <- function(error) sqrt(mean(error^2, na.rm = TRUE))
  RMSE <- Map(function(x) rmse((x - Blue)[]), Filling)
  index <- which.min(unlist(RMSE))
  Filling <- Filling[[index]]
  RMSE <- RMSE[[index]]

  Filling[Filling > 255] <- 255
  Filling[Filling < 0] <- 0

  # fill
  foo <- raster::sampleRegular(Filling, ncell(Filling) * 0.7,
    cells = TRUE
  )
  Filling[foo[, 1]] <- NA

  Blue <- cover(Blue, Filling)

  # Fit trend surface
  r <- .fitTrendSurface(Blue, sampleProportion = 1, np = np)
  GoF <- r[[2]]
  r <- r[[1]]

  if (fact > 1) r <- resample(r, blue)
  compareRaster(r, blue)

  # 10/9/2019
  # A trick to increase the robustness of the fit.
  # It uses a plane to model the DNs near the horizon.
  aux_mask <- mask
  ## filter the estimation
  r[r < 1] <- NA

  ## sample the estimation near the horizon
  aux_mask[] <- raster(EBImage::erode(as.matrix(mask)))[]
  aux_mask <- mask - aux_mask
  aux_r <- r
  aux_r[!aux_mask] <- NA
  ### filter the sample
  thr <- median(aux_r[], na.rm = TRUE) - sd(aux_r[], na.rm = TRUE)
  aux_r[aux_r < thr] <- NA

  ## fit a plane and edit it
  plane <- .fitTrendSurface(aux_r, 1, 1)
  plane <- plane[[1]]
  plane[plane < 1] <- 1
  plane[plane > 254] <- 254
  aux_mask <- z
  aux_mask[!mask] <- NA
  if (getMax(aux_mask) + 10 < 90) {
    aux_mask <- doMask(z, zlim = asAngle(c(
      0,
      getMax(aux_mask) + 10
    )))
  } else {
    aux_mask <- doMask(z)
  }
  plane[aux_mask] <- NA

  ## filter the estimation (again but with other approach)
  r[!mask] <- NA
  thr <- median(r[], na.rm = TRUE) - sd(r[], na.rm = TRUE)
  r[r < thr] <- NA

  ## support the estimation with the plane
  r[!mask] <- NA
  r <- cover(r, plane)

  ## adjust a surface to the supported data
  ### if the result is bad,
  ### it tries with a more rigid model
  m <- doMask(z)
  fun <- function(np) {
    surf <- .fitTrendSurface(r, 0.7, np)
    surf[[1]]
  }
  unlock <- TRUE
  np <- 7
  while (unlock) {
    np <- np - 1
    aux_r <- fun(np)
    r_values <- aux_r[m]
    unlock <- any(r_values < 0)
    if (unlock) unlock <- np > 3
  }

  if (np == 3) {
    thr <- median(r[], na.rm = TRUE)
    aux_r <- fun(6)
    aux_r[aux_r < thr] <- thr
  }

  aux_r[!m] <- NA


  list(skyDN = aux_r, RMSE = RMSE, GoF = GoF)
}
