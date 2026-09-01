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
#' plot(r)
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

  if (any(vals[!is.na(vals)] <= 0)) {
    stop("`flat_field` must contain strictly positive values.",
         call. = FALSE)
  }

  out <- r / flat_field
  out[is.na(out)] <- .get_max(r)
  out[is.na(flat_field)] <- 0

  out
}
