#' Interpolate sky data into a raster and merge it with a sky model raster
#'
#' This function is part of the efforts to automate the method proposed by
#' \insertCite{Lang2010;textual}{rcaiman}. A paper for thoroughly presenting and
#' testing this pipeline is under preparation.
#'
#' @param r [SpatRaster-class]. Typically, the blue channel extracted from a
#'   canopy photograph. The image from which `sky_points` was obtained.
#' @inheritParams ootb_mblt
#' @inheritParams extract_dn
#' @param ootb_sky An object of the class `list` that is the result of calling
#'   [ootb_fit_cie_sky_model()].
#' @param rmax_tune Numeric vector of length one. It must be a positive integer.
#'   It is used to fine tune the `rmax` argument that is computed internally
#'   (see Details).
#'
#' @family Sky Reconstruction Functions
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
#' sky <- ootb_interpolate_and_merge(r, z, a, ootb_sky$sky_points,
#'                              ootb_sky)
#' plot(sky$sky)
#' }
ootb_interpolate_and_merge <- function(r, z, a, sky_points, ootb_sky,
                                  rmax_tune = 1, use_window = TRUE) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  stopifnot(length(rmax_tune) == 1)
  stopifnot(rmax_tune > 0)

  model <- ootb_sky$model
  sky_cie <- sky$sky

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

  i <- grDevices::chull(sky_points$col, nrow(r) - sky_points$row)
  obs <- extract_dn(r, sky_points[i, ])[,3]
  est <- extract_dn(sky_cie, sky_points[i, ])[,3]
  res <- obs - est

  rmax <- pmax(ncol(r)/14,
               ncol(r)/7 * (1 - stats::median(res)/max(obs))) %>% round()
  rmax <- rmax * rmax_tune

  # extract and interpolate
  sky_points <- extract_dn(r, sky_points, use_window = use_window)
  sky <- interpolate_sky_points(sky_points, r, k = k, p = p,
                                rmax = rmax * rmax_tune, col_id = 3)

  # merge
  sky <- sky * (1 - w) + sky_cie * w
  sky <- terra::cover(sky, sky_cie)
  names(sky) <- paste0("Weighted average, ",
                       "w=", round(w, 2), ", ",
                       "k=", k, ", ",
                       "p=", round(p, 2), ", ",
                       "rmax=", rmax)

  list(sky = sky,
       w = w,
       k=  k,
       p = p,
       rmax = rmax)
}

