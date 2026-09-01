#' Compute mean
#'
#' @description
#' Compute the mean value of raster pixels within each segment. This is a
#' lightweight alternative to [compute_features()] when only segment means are
#' required.
#'
#' @param r numeric [terra::SpatRaster-class] with one layer. Pixel values used
#'   to compute object mean.
#' @param seg single-layer [terra::SpatRaster-class]. Segmentation map
#'   of `r`, typically created with functions such as [healpix_segmentation()],
#'   [equal_area_segmentation()] or [sky_grid_segmentation()], but any raster with
#'   integer segment labels is accepted.
#' @param return character of length one. Either `"raster"` (default) or
#'   `"vector"`, controlling whether to return a map with per-segment values or
#'   a named vector (one value per segment).
#' @param ignore_label_0 logical of length one. If `TRUE`, the segment labeled
#'   `0` is ignored.
#'
#' @return If `return = "raster"`, a [terra::SpatRaster-class] where each pixel
#'   holds its segment’s feature value. If `return = "vector"`, a named numeric
#'   (or logical) vector with one value per segment.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' print(compute_mean(r$Cyan, g, return = "vector"))
#' plot(compute_mean(r$Cyan, g, return = "raster"))
#' }
compute_mean <- function(r, seg,
                         return = "raster",
                         ignore_label_0 = TRUE) {

  .assert_single_layer(r, allow_na = TRUE)
  .assert_single_layer(seg)
  .assert_same_geom(r, seg)

  .assert_choice(return, c("vector", "raster"))
  .check_vector(ignore_label_0, "logical", 1)

  out <- compute_mean_cpp(
    terra::values(r),
    terra::values(seg),
    ignore_label_0
  )

  if (return == "vector") {
    id <- out$id
    out <- out$mean
    names(out) <- id
    return(out)
  } else {
    feature_image(seg, as.data.frame(out), "mean")
  }
}
