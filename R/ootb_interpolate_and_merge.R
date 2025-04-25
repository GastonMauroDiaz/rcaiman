#' Interpolate sky data into a raster and merge it with a sky model raster
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly presenting and
#' testing this pipeline is under preparation.
#'
#' @param r [SpatRaster-class]. Typically, the blue channel extracted from a
#'   canopy photograph. The image from which `sky_points` was obtained.
#' @inheritParams sky_grid_segmentation
#' @inheritParams extract_dn
#' @param ootb_sky An object of the class `list` that is the result of calling
#'   [ootb_fit_cie_sky_model()].
#' @param size Numeric vector of length one.
#'
#' @return An object of class [SpatRaster-class].
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path), z, a)
#'
#' sky <- ootb_interpolate_and_merge(r, z, a, ootb_sky$sky_points, ootb_sky)
#'
#' plot(sky$sky)
#'
#' # See fit_cie_sky_model() for details on below file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- rbind(sky_points, ootb_sky$sky_points[, c("row", "col")])
#'
#' sky <- ootb_interpolate_and_merge(r, z, a, sky_points, ootb_sky)
#'
#' plot(sky$sky)
#' }
ootb_interpolate_and_merge <- function(r, z, a, sky_points, ootb_sky,
                                  size = 100,
                                  use_window = TRUE) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  model <- ootb_sky$model
  sky_cie <- ootb_sky$sky

  # read metrics
  r2 <- ootb_sky$model_validation$r_squared
  error <- stats::median(abs(1 - ootb_sky$model_validation$obs /
                               ootb_sky$model_validation$pred))

  # obtain parameters
  w <-  1 - stats::plogis(error, 0.042, 0.05) + stats::plogis(0, 0.042, 0.05)
  if (!is.numeric(r2)) r2 <- 0
  k <- 2 + round(r2 * 10)
  p <- abs(model$coef[1]) + log(model$coef[3]+1)
  p <- 2 + sqrt(p)

  stopifnot(k <= nrow(sky_points))

  # L2010
  sky_points2 <- expand_sky_points(r, z, a, sky_points,
                                   sky_model =  sky_cie,
                                   k = k,
                                   p = p,
                                   w = w,
                                   rule = "any",
                                   chi_max = 20,
                                   size = size)
  sky <-  interpolate_sky_points(sky_points2[!sky_points2$initial, ], r,
                                 k = 1, p = 1,
                                 rmax = ncol(r)/size + 1,
                                 col_id = "dn")
  sky[is.na(z)] <- NA

  names(sky) <- "Sky vault"

  list(sky = sky,
       w = w,
       k=  k,
       p = p)
}

