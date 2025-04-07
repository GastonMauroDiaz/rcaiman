#' Expand non-circular
#'
#' Expand a non-circular hemispherical photograph.
#'
#' @param caim [SpatRaster-class]. The output of [read_caim()] or
#'   [read_caim_raw()].
#' @inheritParams ootb_mblt
#' @param zenith_colrow Numeric vector of length two. Raster coordinates of the
#'   zenith. See [calc_zenith_colrow()].
#'
#' @family Lens Functions
#'
#' @return An object of class [SpatRaster-class] that is the result of adding
#'   margins (`NA` pixels) to the argument `caim`. The zenith point depicted in
#'   the picture should be in the center of the image or very close to it.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#'  # ====================================================================
#'  # Non-circular Fisheye images from a Smartphone with an Auxiliary Lens
#'  # (Also applicable to Non-circular images from DSLR cameras)
#'  # ====================================================================
#'
#'  path <- system.file("external/APC_0581.jpg", package = "rcaiman")
#'  caim <- read_caim(path)
#'  z <- zenith_image(2132/2,  c(0.7836, 0.1512, -0.1558))
#'  a <- azimuth_image(z)
#'  zenith_colrow <- c(1063, 771)/2
#'  caim <- expand_noncircular(caim, z, zenith_colrow)
#'  plot(caim$Blue, col = seq(0,1,1/255) %>% grey())
#'  m <- !is.na(caim$Red) & !is.na(z)
#'  plot(m, add = TRUE, alpha = 0.3, legend = FALSE)
#'
#'
#'  # ============================
#'  # Restricted View Canopy Photo
#'  # ============================
#'
#'  path <- system.file("external/APC_0020.jpg", package = "rcaiman")
#'  caim <- read_caim(path)
#'  plot(caim)
#'  caim <- normalize_minmax(caim)
#'  diameter <- calc_diameter(lens(), sqrt(nrow(caim)^2 + ncol(caim)^2)/2, 90)
#'  z <- zenith_image(diameter, lens())
#'  caim <- expand_noncircular(caim, z, c(ncol(caim)/2, nrow(caim)/2))
#'  m <- !is.na(caim$Red)
#'  a <- azimuth_image(z)
#'  caim[!m] <- 0
#'  z <- normalize_minmax(z, 0, 90) * 30.15 #60.3º diagonal FOV according to metadata
#'  plot(caim$Blue, col = seq(0,1,1/255) %>% grey())
#'  m <- !is.na(caim$Red) & !is.na(z)
#'  plot(m, add = TRUE, alpha = 0.3, legend = FALSE)
#' }
expand_noncircular <-  function (caim, z, zenith_colrow) {
  .is_single_layer_raster(z, "z")
  stopifnot(class(zenith_colrow) == "numeric")
  stopifnot(length(zenith_colrow) == 2)

  zenith_xy <- c(zenith_colrow[1], nrow(caim) - zenith_colrow[2])
  delta_x <-  zenith_xy[1] - (ncol(caim) / 2)
  delta_y <-  zenith_xy[2] - (nrow(caim) / 2)
  #In which quadrant is the zenith?
  # (-)(+)|(+)(+)
  #----------------
  # (-)(-)|(+)(-)
  #
  center <- ncol(z) / 2
  xmn <- center - ((ncol(caim)/2) + delta_x)
  xmx <- center + ((ncol(caim)/2) - delta_x)
  ymn <- center - ((nrow(caim)/2) + delta_y)
  ymx <- center + ((nrow(caim)/2) - delta_y)
  e <- terra::ext(xmn, xmx, ymn, ymx)
  r <- terra::deepcopy(caim)
  terra::ext(r) <- e

  r <- terra::extend(r, z)
  terra::ext(r) <- terra::align(terra::ext(r), z)
  r <- terra::resample(r, z)
  terra::ext(r) <- terra::ext(0, ncol(r), 0, nrow(r))
  r
}
