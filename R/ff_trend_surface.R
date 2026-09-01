#' Generate flat-field using a trend surface model
#'
#' @description
#' Generate a flat-field image from dense angular radiance samples using
#' a polynomial trend surface model.
#'
#' @details
#' Radiance values are selected for a given `aperture` and `channel`, mapped from
#' angular sampling positions (`zenith`, `azimuth`, in degrees) to raster
#' coordinates, and fitted with [fit_trend_surface()] `np = 6`.
#'
#' The resulting surface represents the relative spatial sensitivity of the
#' imaging system under uniform radiance conditions.
#'
#' @note
#' A warning is issue if the coefficient of determination of the fit is lower
#' than 0.85.
#'
#' @param data named list. List names must be coercible to numeric values
#'   representing f-numbers. Each aperture must have dense flat-field data apt
#'   to build a flat-field image. The expected structure of each component is
#'   another list with components:
#'   \itemize{
#'     \item `direction`: data frame with numeric columns `zen` and
#'     `azi` (radian), describing the direction as in [add_directions()].
#'     \item `relative_radiance`: data frame with at least one numeric column
#'     with the same number of rows as `direction`. Columns correspond to CFA
#'     channels as it is conceptualized in [split_cfa()]).
#'   }
#'   Angular values must be finite. Radiance values must be finite and
#'   non-negative.
#' @param aperture character vector of length one. Aperture value available in
#'   `data`.
#' @param channel character vector. Name of colum or columns of the
#'   `relative_radiance` data frame in `data`.
#'
#' @inheritParams add_directions
#' @inheritParams fit_trend_surface
#'
#'
#' @return
#'   Numeric [terra::SpatRaster-class] containing the flat-field surface or
#'   surfaces named as `channel`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw2.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#'
#' rspec <- get_radiometry_spec(spec, "flat_field_reconstruction")
#' z <- zenith_image(spec$horizon_radius * 2, spec$parameters)
#' a <- azimuth_image(z)
#' ff <- ff_trend_surface(z, a, rspec$data, "5", "Cyan")
#' plot(ff)
#' }
ff_trend_surface <- function(z, a, data, aperture, channel) {
  .validate_flat_field_data(data)
  .check_vector(aperture, "character", 1)
  .check_vector(channel, "character")

  available_apertures <- names(data$relative_radiance)

  if (!aperture %in% available_apertures) {
    stop(
      sprintf(
        "Invalid `aperture` ('%s'). Available apertures: %s",
        aperture,
        paste(available_apertures, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  available_channels <- names(data$relative_radiance[[aperture]])


  row_col <- row_col_from_zen_azi(z, a, data$direction$zen, data$direction$azi)

  l <- list()
  for (u in channel) {

    if (!u %in% available_channels) {
      stop(
        sprintf(
          "Invalid `channel` ('%s') for aperture '%s'. Available channels: %s",
          u,
          aperture,
          paste(available_channels, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    sp <- data.frame(row_col, data$relative_radiance[[aperture]][[u]])
    names(sp)[3] <- u
    ff <- fit_trend_surface(sp, col_id = 3, z, np = 6, TRUE, return = "all")
    if (ff$r2 < 0.85)
      warning(
         sprintf("Coefficient of determination is equal to %s", ff$r2)
        )
    ff <- ff$raster
    ff[is.na(z)] <- NA
    names(ff) <- u
    l <- c(l, ff)
  }
  rast(l)
}
