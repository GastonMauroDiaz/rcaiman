#' Calc sky-normalized green difference (SNGD)
#'
#' @inheritParams membership_to_color
#' @inheritParams ootb_mblt
#' @param bin [SpatRaster-class]. It should be the result of binarizing an
#'   enhaced canopy image and post-processing it to remove microgaps. See
#'   the example.
#' @param angle Numeric vector of length one. Zenith angle in degrees. The
#'   canopy below this value of zenith angle will not be analyzed.
#' @param kern_size Numeric vector of length one. Controls the morphological
#'   operation used to remove pixels on the canopy silhouette contour, which are likely
#'   to be mixed pixeles.
#'
#' @family Tool Functions
#' @returns Numeric vector of length one.
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- .caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' r <- caim$Blue
#'
#' bin <- regional_thresholding(r, rings_segmentation(z, 30),
#'                              method = "thr_isodata")
#' g <- sky_grid_segmentation(z, a, 10)
#' sun_coord <- extract_sun_coord(r, z, a, bin, g)
#' sun_coord$zenith_azimuth
#'
#' .a <- azimuth_image(z, orientation = sun_coord$zenith_azimuth[2]+90)
#' seg <- sectors_segmentation(.a, 180) * rings_segmentation(z, 30)
#' bin <- regional_thresholding(r, seg, method = "thr_isodata")
#'
#' mx <- optim_normalize(caim, bin)
#' caim <- normalize(caim, mx = mx, force_range = TRUE)
#' ecaim <- enhance_caim(caim, m, polarLAB(50, 17, 293))
#' bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))
#' bin <- apply_thr(terra::focal(bin, 3, median), 0.5)
#'
#' calc_sngd(.caim, z, a, bin, kern_size = 3)
#' }
calc_sngd <- function(caim, z, a, bin,
                                angle = 60,
                                kern_size = 3) {
  stopifnot(!.is_even(kern_size))
  stopifnot(.is_whole(kern_size))
  .this_requires_EBImage()
  terra::compareGeom(caim, z)
  terra::compareGeom(caim, bin)
  .is_logic_and_NA_free(bin, "bin")

  # remove mixed pixels
  kern <- EBImage::makeBrush(kern_size, "box")
  bin <- EBImage::dilate(as.array(bin), kern) %>%terra::setValues(bin, .)

  # mask plant pixels that are lighter than the others plant pixels
  bin2 <- apply_thr(caim$Green, thr_isodata(caim$Green[m]))
  bin2[bin] <- 0

  # find clusters of light plant pixel considering solid angle
  seg <- polar_qtree(bin2, z, a, scale_parameter = 0)
  bin2 <- extract_feature(bin2, seg, mean) > 0.9

  # clip sky pixels and every pixels outside the AOI
  m <- !bin & select_sky_vault_region(z, 0, angle)
  bin2[!m] <- NA

  if (any(as.logical(bin2[m]))) {
    v <- tapply(caim$Green[!is.na(bin2)], bin2[!is.na(bin2)], stats::median)
    mx <- max(caim$Blue[!is.na(z)])
    sunlit_canopy_index <- unname((v[2] - v[1])/mx * 100)
  } else {
    sunlit_canopy_index <- 0
  }

  sunlit_canopy_index
}
