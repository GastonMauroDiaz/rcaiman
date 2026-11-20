#' Apply flat-field correction
#'
#' @description
#' Correct an image for spatially varying system response using a flat-field
#' image.
#'
#' @details
#' Flat-field correction compensates for spatially varying system response
#' (e.g., vignetting) by dividing the observed radiance image by a flat-field
#' image representing the relative spatial sensitivity of the imaging system
#' under uniform radiance conditions.
#'
#' If the response were perfectly uniform, the flat-field image would be
#' constant and equal to one.
#'
#' The correction is defined as:
#'
#' \deqn{r_{corrected} = r / F}
#'
#' where \eqn{F} is the flat-field image representing the relative spatial
#' response of the imaging system.
#'
#' If multiple bands are present, correction is applied band-wise. The flat
#' field must contain strictly positive and finite values.
#'
#' @param r [terra::SpatRaster-class]
#' @param flat_field [terra::SpatRaster-class]. Flat-field image with the same
#'   spatial dimensions and number of bands as `r`. Must contain strictly
#'   positive values or `NA`.
#'
#' @return Object of the same class as `r`, with values divided by
#'   `flat_field`.
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(30, 20),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#'  )
#' g <- set_internal_geometry(c(60, 40), geometry)
#' r <- g$off_axis
#' ff <- g$off_axis
#' corrected_r <- apply_flat_field(r, ff)
#' plot(corrected_r)
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/DSCN4248.tif", package = "rcaiman")
#' metadata <- read.delim(gsub("tif", "tsv", path))
#' sensor_data <- read_caim(path)
#'
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#'
#' my_entry <- read_registry_snapshot(path)
#' raw <- get_signature(my_entry, "raw")
#' cygm <- get_band_spec(raw, "cygm")
#'
#' caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
#' caim
#' r <- caim$Cyan
#'
#' raw_geometry <- get_geometry_spec(cygm, "simple_method")
#' g <- set_internal_geometry(raw$dim, raw_geometry)
#'
#' raw_radiometry <- get_radiometry_spec(raw_geometry, "photometric_sphere")
#'
#' aperture <- as.character(metadata$aperture)
#' ff <- trend_surface_image(r, raw_radiometry$parameters[[aperture]]$Cyan,
#'                           m = !is.na(g$off_axis))
#' corrected_r <- apply_flat_field(r, ff)
#' plot(corrected_r)
#' }
apply_flat_field <- function(r, flat_field) {

  .assert_spatraster(r)
  .assert_spatraster(flat_field)
  .assert_same_geom(r, flat_field)

  if (terra::nlyr(r) != terra::nlyr(flat_field)) {
    stop("`r` and `flat_field` must have the same number of layers.",
         call. = FALSE)
  }

  vals <- terra::values(flat_field, mat = FALSE)

  if (any(!is.finite(vals[!is.na(vals)]) | vals[!is.na(vals)] <= 0)) {
    stop("`flat_field` must contain finite, strictly positive values.",
         call. = FALSE)
  }

  out <- r / flat_field
  out[is.na(out)] <- .get_max(r)
  out[is.na(flat_field)] <- 0

  out
}
