#' Build internal viewing directions
#'
#' @description
#' Describe each pixel's viewing direction relative to the optical axis and
#' the top of the sensor. The acquisition orientation is not involved.
#' The radial model interface is stable; the anisotropic model is experimental.
#'
#' @param dim numeric vector of length two with positive integer-like values,
#'   `c(width, height)` in pixels.
#' @param geometry named list describing the internal geometry. No particular
#'   class is required; additional fields are ignored. The fields used are:
#'   \describe{
#'     \item{`model`}{either `"radial_projection"` or
#'       `"anisotropic_projection"`.}
#'     \item{`center_of_vision`}{integer-like numeric vector of length two,
#'       `c(x, y_down)`, locating the optical axis. Coordinates start at the
#'       upper-left image corner and increase rightward and downward. This is
#'       the **snapped** center: integer coordinates are cell boundaries; the
#'       first pixel center is `c(0.5, 0.5)`. See [y_down_from_y()].}
#'     \item{`parameters`}{model coefficients. For the radial model, a nonempty
#'       finite numeric vector (`projection_coef`), in increasing order of
#'       degree and without an intercept. For the anisotropic model, a named
#'       list described below.}
#'     \item{`normalization_radius`}{positive numeric value. Radius in pixels
#'       corresponding to `normalization_angle`. Required for the radial model.}
#'     \item{`normalization_angle`}{positive off-axis angle in degrees, at most
#'       180, where the radial polynomial must equal one (tolerance 0.01).
#'       Underestimation beyond floating-point roundoff is rejected; use
#'       [fix_projection_coef()] explicitly. See [validate_projection_function()].
#'       Required for the radial model. Declaring this reference does not
#'       rescale the coefficients.}
#'     \item{`max_off_axis_angle`}{positive angular limit in degrees, at most
#'       180, for radial inversion. Defaults to 90 when absent. Pixels beyond
#'       its projected radius have `NA` in both layers. This is a calculation
#'       limit, independent of normalization, not a claim about calibration
#'       coverage.}
#'   }
#'   The anisotropic model ignores the normalization fields and angular limit.
#'   Its `parameters` list must contain exactly:
#'   \describe{
#'     \item{`off_axis_coef`}{nonempty finite numeric vector of coefficients
#'       predicting off-axis angle from the transformed radius, in increasing
#'       order of degree and without an intercept.}
#'     \item{`sin_internal_azimuth_coef`}{finite numeric vector of length two
#'       predicting the sine of internal azimuth from `sin(phi)` and
#'       `cos(phi)`, in that order and without an intercept.}
#'     \item{`cos_internal_azimuth_coef`}{finite numeric vector of length two
#'       predicting the cosine of internal azimuth from the same predictors.}
#'     \item{`atan_scale`}{positive numeric value. Hyperparameter controlling
#'       the transformation of relative radius.}
#'   }
#'   See *Details* for the equations and angular conventions of both models.
#'
#' @details
#' A viewing ray represents a line of sight from the entrance pupil into the
#' scene. Internal geometry describes its direction relative to the optical
#' axis and the top of the sensor.
#'
#' The off-axis angle is measured from the optical axis, just as off-nadir
#' and off-zenith angles are measured from nadir and zenith, respectively.
#' When the optical axis points vertically upward, the off-axis angle equals
#' the zenith angle.
#'
#' Internal azimuth is measured clockwise from the top of the sensor in its
#' reference landscape disposition: top = 0, right = 90, bottom = 180, and
#' left = 270 deg. Rotating the device does not change these internal
#' directions. The output grid retains the dimensions supplied in `dim`.
#'
#' Both models use offsets `dx, dy_down` from `center_of_vision`, measured
#' rightward and downward at pixel centers. Their distance to the center is
#' \eqn{r = \sqrt{dx^2 + dy_{down}^2}} in pixels. The equations below use
#' radians; the returned layers use degrees.
#'
#' ### Radial projection
#'
#' `r lifecycle::badge("stable")`
#'
#' This model assumes circular symmetry: equal distances from the center
#' correspond to equal off-axis angles. The projection polynomial predicts
#' relative radius from off-axis angle:
#' \deqn{\rho = \frac{r}{r_{norm}} = f(\theta') =
#'   \sum_{k=1}^{p} a_k (\theta')^k}
#' Here, \eqn{r_{norm}} is `normalization_radius`, and the coefficients
#' \eqn{a_1, \ldots, a_p} are the `projection_coef` vector supplied through
#' `parameters`. Its length sets the polynomial degree. The normalization
#' reference satisfies \eqn{f(\theta'_{norm}) = 1}, where
#' \eqn{\theta'_{norm}} is `normalization_angle` expressed in radians.
#' See [normalize_projection_coef()] to change this reference and the
#' corresponding radius in pixels, including for the standard projections.
#'
#' To recover off-axis angle at a pixel, invert this polynomial on
#' `[0, max_off_axis_angle]` (degrees), where it must be strictly increasing.
#' [calc_relative_radius()] and [calc_off_axis_angle()] evaluate the same
#' model in either direction for numeric vectors.
#' Relative radii above one are allowed within this interval. Internal azimuth
#' follows directly from the image offsets:
#' \deqn{\alpha = \operatorname{atan2}(dx, -dy_{down}) \bmod 2\pi}
#'
#' ### Anisotropic projection
#'
#' `r lifecycle::badge("experimental")`
#'
#' This model predicts off-axis angle from a transformed relative radius and
#' internal azimuth from the image-plane angle. Normalize distance by the
#' image semidiagonal:
#' \deqn{\rho = \frac{r}{\sqrt{(w/2)^2 + (h/2)^2}},
#'   \qquad t = \arctan(\rho / s)}
#' Here, \eqn{w} and \eqn{h} are the width and height in `dim`, and
#' \eqn{s} is the hyperparameter `atan_scale`. The off-axis polynomial is:
#' \deqn{\theta' = \sum_{k=1}^{p} b_k t^k}
#' The coefficients \eqn{b_1, \ldots, b_p} are `off_axis_coef`.
#'
#' The image-plane angle \eqn{\phi = \operatorname{atan2}(dy_{down}, dx)}
#' starts at the right and increases clockwise. Fit the sine and cosine of
#' internal azimuth to avoid the discontinuity where angles wrap around:
#' \deqn{S = u_1\sin(\phi) + u_2\cos(\phi)}
#' \deqn{C = v_1\sin(\phi) + v_2\cos(\phi)}
#' The vectors \eqn{(u_1, u_2)} and \eqn{(v_1, v_2)} are
#' `sin_internal_azimuth_coef` and `cos_internal_azimuth_coef`, respectively.
#' Recover internal azimuth, measured clockwise from the top, with:
#' \deqn{\alpha = \operatorname{atan2}(S, C) \bmod 2\pi}
#'
#' @section Masks in the internal frame:
#' Compare the angular layers directly to select regions tied to the optical
#' system, such as a sector obstructed by a camera mount. These masks use
#' degrees and do not depend on acquisition orientation. For sky regions
#' relative to zenith and north, see [direction_image()].
#'
#' Include `!is.na()` for each angle used so that missing directions receive
#' `FALSE`. Use `&` for the intersection of conditions and `|` for their union.
#' A sector crossing zero requires two intervals. For example, to exclude an
#' illustrative obstruction beyond 70 deg off axis and within 15 deg of the
#' top of the sensor:
#'
#' ```r
#' a <- g$internal_azimuth
#' valid <- !is.na(g$off_axis) & !is.na(a)
#' obstruction <- valid & g$off_axis >= 70 & (a >= 345 | a <= 15)
#' m <- valid & !obstruction
#' ```
#'
#' Explicit comparisons determine whether limits are included (`>=`, `<=`)
#' or excluded (`>`, `<`). A selection may legitimately include all available
#' pixels or none. Use the resulting logical mask with [rem_unselected()],
#' [paint_with_mask()], or other functions accepting `m`.
#'
#' @note
#' For the anisotropic model, The declared center and image semidiagonal must
#' match those used to fit the model and off-axis predictions are not clipped to a
#' hemisphere.
#'
#' @return numeric [terra::SpatRaster-class] with two layers on the pixel grid
#'   established by [set_rcaiman_geometry()]:
#'   \describe{
#'     \item{`off_axis`}{angle from the optical axis, in degrees.}
#'     \item{`internal_azimuth`}{clockwise angle from the top of the sensor,
#'       in degrees, in `[0, 360)`.}
#'   }
#'
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(40, 30),
#'   parameters = c(0.6380, 0.0306, -0.0199),
#'   normalization_radius = 25,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' plot(g)
#'
#' # Normalize at 60 degrees while calculating directions up to 90 degrees.
#' geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(40, 30),
#'   parameters = 1 / (pi / 3),
#'   normalization_radius = 20,
#'   normalization_angle = 60
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#'
#' a <- g$internal_azimuth
#' valid <- !is.na(g$off_axis) & !is.na(a)
#' obstruction <- valid & g$off_axis >= 70 & (a >= 345 | a <= 15)
#' m <- valid & !obstruction
#' plot(m)
#'
#' # Illustrative internal coefficients, not an instrument calibration.
#' geometry <- list(
#'   model = "anisotropic_projection",
#'   center_of_vision = c(40, 30),
#'   parameters = list(
#'     off_axis_coef = 1,
#'     sin_internal_azimuth_coef = c(0, 1),
#'     cos_internal_azimuth_coef = c(-1, 0),
#'     atan_scale = 0.6
#'   )
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
set_internal_geometry <- function(dim, geometry) {
  .check_vector(dim, "integerish", 2, sign = "positive")
  .check_internal_geometry(geometry)
  center_of_vision <- geometry[["center_of_vision"]]
  model <- geometry[["model"]]
  parameters <- geometry[["parameters"]]
  normalization_radius <- geometry[["normalization_radius"]]
  max_off_axis_angle <- if ("max_off_axis_angle" %in% names(geometry))
    geometry[["max_off_axis_angle"]] else 90

  g <- terra::rast(ncols = dim[1], nrows = dim[2], nlyrs = 2)
  g <- set_rcaiman_geometry(g)
  xy <- terra::xyFromCell(g, seq_len(terra::ncell(g)))
  dx <- xy[, 1] - center_of_vision[1]
  dy <- dim[2] - xy[, 2] - center_of_vision[2]

  angles <- switch(model,
    radial_projection = .internal_radial_directions(
      dx, dy, parameters, normalization_radius, max_off_axis_angle
    ),
    anisotropic_projection = .internal_anisotropic_directions(
      dx, dy, dim, parameters
    )
  )
  terra::values(g) <- angles * 180 / pi
  names(g) <- c("off_axis", "internal_azimuth")
  # Rounding immediately below 2*pi can reach 360 during unit conversion.
  azi <- terra::values(g[[2]], mat = FALSE)
  azi[!is.na(azi) & azi >= 360] <- 0
  g[[2]] <- terra::setValues(g[[2]], azi)
  g
}

.internal_radial_directions <- function(dx, dy, parameters, radius, max_angle) {
  relative_radius <- sqrt(dx^2 + dy^2) / radius
  off_axis <- .calc_off_axis_angle(relative_radius, parameters, max_angle)
  azimuth <- atan2(dx, -dy) %% (2 * pi)
  azimuth[is.na(off_axis) | relative_radius == 0] <- NA_real_
  cbind(off_axis, azimuth)
}

.internal_anisotropic_directions <- function(dx, dy, dim, parameters) {
  polar <- .image_polar_predictors(dx, dy, dim[1], dim[2])
  atan_radius <- atan(polar$relative_radius / parameters$atan_scale)
  off_axis <- numeric(length(atan_radius))
  for (degree in seq_along(parameters$off_axis_coef)) {
    off_axis <- off_axis + parameters$off_axis_coef[degree] * atan_radius^degree
  }
  phi <- polar$image_angle
  s <- parameters$sin_internal_azimuth_coef[1] * sin(phi) +
    parameters$sin_internal_azimuth_coef[2] * cos(phi)
  c <- parameters$cos_internal_azimuth_coef[1] * sin(phi) +
    parameters$cos_internal_azimuth_coef[2] * cos(phi)
  azimuth <- atan2(s, c) %% (2 * pi)
  azimuth[!is.na(azimuth) & azimuth >= 2 * pi] <- 0
  azimuth[is.na(dx) | is.na(dy) | (dx == 0 & dy == 0) | (s == 0 & c == 0)] <- NA_real_
  cbind(off_axis, azimuth)
}

#' Orient viewing directions
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Express internal viewing directions as zenith and geographic azimuth angles
#' for the orientation of an acquisition.
#'
#' @details
#' Start with the device level, looking south, with its sensor in landscape
#' disposition with its top upward. Describe the orientation as three movements
#' of device mounted on a tripod:
#' \enumerate{
#'   \item Rotate the sensor clockwise around the optical axis, viewed from
#'     behind the device toward the scene (`sensor_rotation`). For example,
#'     rotations of 90 or 270 deg specify the two portrait dispositions.
#'   \item Turn horizontally to the right from south (`azimuth_rotation`): 0 deg
#'     faces south, 90 deg west, 180 deg north, and 270 deg east.
#'   \item Raise the optical axis above the horizontal (`elevation`).
#'     Its elevation is independent of the sensor's landscape or portrait
#'     disposition.
#' }
#' This sequence defines the orientation; it does not prescribe the order
#' in which the operator adjusts the device in the field.
#'
#' The defaults describe a vertical upward acquisition with the top of the
#' sensor toward north. In this configuration, zenith angle equals off-axis
#' angle, and geographic azimuth is the negative of internal azimuth, modulo
#' 360 deg. North is at the top, east at the left, and west at the right.
#' For a hinge acquisition at 57.5 deg from zenith, use `elevation = 32.5`.
#'
#' `azimuth_rotation` describes the horizontal turn from south. When the
#' optical axis is vertical, it still determines the orientation of the
#' sensor, although the axis itself has no geographic azimuth.
#'
#' Acquisition orientation allows the same viewing directions to be expressed
#' relative to the local vertical and north. These external angles describe
#' where the instrument is looking, without specifying the position from which
#' the image was acquired. Images acquired at different positions can therefore
#' share the same external directions while sampling different paths through
#' the canopy.
#'
#' Directions are transformed using unit vectors, as in
#' [external_directions_from_internal()]. The output contains only the external
#' directions, on the same grid. The internal geometry remains unchanged and
#' can be reused with another acquisition orientation. Directions below the
#' horizon have zenith angles greater than 90 deg and are retained.
#'
#' @section Masks in the external frame:
#' Compare `d$zenith` and `d$azimuth` to select sky regions independently of
#' acquisition orientation. The layers use degrees. For example, a zenith
#' interval with both boundaries included is:
#'
#' ```r
#' z <- d$zenith
#' m <- !is.na(z) & z >= 20 & z <= 70
#' ```
#'
#' For a geographic sector crossing north, combine two azimuth intervals:
#'
#' ```r
#' a <- d$azimuth
#' north <- !is.na(a) & (a >= 345 | a <= 15)
#' m <- m & north
#' ```
#'
#' These comparisons produce logical masks without `NA`. A zenith-only mask
#' does not require an azimuth, so it can include the vertical direction.
#' A mask requiring azimuth excludes pixels where azimuth is undefined.
#' To exclude a region, combine its complement with an availability mask,
#' for example `!is.na(z) & !m`; using `!m` alone would also select pixels
#' without directions. For regions tied to the equipment, use internal angles
#' as illustrated in [set_internal_geometry()].
#'
#' @section Describing an imperfectly levelled acquisition:
#' In upward-looking photography, the usual procedure is to mount the device on
#' a tripod, orient the lens upward, and level it so that the optical axis is
#' vertical. A bubble level indicates both how far the device is tilted and
#' which way to correct it, even though we rarely express the latter as an angle.
#'
#' When describing an imperfectly levelled acquisition, both pieces of
#' information are needed: the inclination from vertical and its direction.
#' The orientation of the sensor must also be specified.
#'
#' The arguments of `direction_image()` describe these properties through a
#' sequence that starts with the device horizontal and looking south. This is
#' convenient for inclined acquisitions, such as photography at the hinge
#' angle. For nearly vertical acquisitions, the same orientation can be easier
#' to describe starting from an upward-looking position.
#'
#' Start with the optical axis pointing vertically upward. Specify the
#' direction toward which the top of the sensor faces---the side where the
#' flash is located on a typical SLR camera. With `sensor_rotation = 0`,
#' `azimuth_rotation = 0` places the top toward north, and 90 deg toward east:
#'
#' ```r
#' sensor_top_azimuth <- 0
#' d <- direction_image(g, azimuth_rotation = sensor_top_azimuth)
#' ```
#'
#' Then describe a departure from that position using an inclination and its
#' direction:
#'
#' ```r
#' inclination <- 3
#' inclination_direction <- 90
#'
#' d <- direction_image(
#'   g,
#'   sensor_rotation = inclination_direction,
#'   azimuth_rotation = sensor_top_azimuth + inclination_direction,
#'   elevation = 90 - inclination
#' )
#' ```
#'
#' Here, `inclination` is the angular departure from vertical, in degrees.
#' `inclination_direction` indicates which way the optical axis leans relative
#' to the initial upward-looking image: 0 degrees toward its bottom, 90 toward
#' its right, 180 toward its top, and 270 toward its left. Thus, this example
#' describes a 3-degree inclination toward the right of that initial image.
#'
#' The two rotation arguments work together. With zero inclination, their
#' contributions from `inclination_direction` cancel, leaving the original
#' upward-looking orientation. With nonzero inclination, their combination
#' selects the direction of the lean. `sensor_rotation` therefore has a role
#' beyond specifying portrait or landscape disposition.
#'
#' These variables describe an orientation; they do not estimate it from the
#' image. The inclination and its direction must come from measurements or an
#' alignment procedure.
#'
#' When aligning acquisitions by reprojection, describe each acquisition with
#' its own `d`, then use the same destination geometry, dimensions, and
#' orientation in [reproject_image()]. This keeps the orientation adjustment
#' associated with each acquisition while ensuring that corresponding output
#' cells represent the same viewing directions.
#'
#' [fit_image_orientation()] can provide an initial orientation from
#' corresponding landmarks. To assess and refine it, reproject both images to
#' the common destination, load them as overlaid layers in
#' [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects),
#' and use flicker to compare their alignment. Adjust the orientation based on
#' that comparison, regenerate `d`, and reproject the original image again.
#' Repeat this visual comparison as needed, keeping the destination fixed.
#'
#' @param g numeric [terra::SpatRaster-class] with two layers
#'   named `off_axis` and `internal_azimuth`, in that order, as returned by
#'   [set_internal_geometry()].
#' @param sensor_rotation numeric vector of length one. Clockwise rotation of
#'   the sensor in degrees. Zero represents landscape; 90 and 270 deg are the
#'   two portrait dispositions. Negative angles and angles greater than 360 deg
#'   are allowed. Full turns do not change the resulting orientation.
#' @param azimuth_rotation numeric vector of length one. Horizontal turn to the
#'   right from south, in degrees. Negative angles and angles greater than 360
#'   deg are allowed. Full turns do not change the resulting orientation.
#' @param elevation numeric vector of length one. Elevation of the optical
#'   axis above the horizontal in degrees, from -90 to 90 deg. Positive values
#'   point upward; 90 deg points to zenith and -90 deg to nadir.
#'
#' @return numeric [terra::SpatRaster-class] with the same grid as
#'   `g` and two layers:
#'   \describe{
#'     \item{`zenith`}{angle from the upward vertical, in degrees, in
#'       `[0, 180]`.}
#'     \item{`azimuth`}{geographic azimuth, in degrees, in `[0, 360)`:
#'       north = 0, east = 90, south = 180, and west = 270 deg. At zenith
#'       and nadir, azimuth is `NA`.}
#'   }
#'   Missing internal directions produce missing external directions.
#'
#' @note
#' A pixel pointing to zenith or nadir has a defined zenith angle (0 or
#' 180 deg), but no defined azimuth. Its `azimuth` value is therefore `NA`,
#' even within a region of valid observations. This `NA` describes an
#' undefined angular coordinate, not missing image data. Check the `zenith`
#' layer to distinguish this case from a missing direction.
#'
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(40, 30),
#'   parameters = 1 / (pi / 3),
#'   normalization_radius = 20,
#'   normalization_angle = 60
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d_upward <- direction_image(g)
#' d_hinge <- direction_image(g, elevation = 32.5)
#' d_hinge_west <- direction_image(g, azimuth_rotation = 90, elevation = 32.5)
#' d_portrait <- direction_image(g, sensor_rotation = 90, elevation = 32.5)
#' plot(d_hinge)
#'
#' z <- d_hinge$zenith
#' a <- d_hinge$azimuth
#' m <- !is.na(z) & z >= 20 & z <= 70
#' north <- !is.na(a) & (a >= 345 | a <= 15)
#' plot(m & north)
direction_image <- function(g,
                            sensor_rotation = 0,
                            azimuth_rotation = 0,
                            elevation = 90) {
  .assert_spatraster(g)
  if (!identical(names(g),
                 c("off_axis", "internal_azimuth"))) {
    stop("`g` must have two layers named ",
         "`off_axis` and `internal_azimuth`, in that order.", call. = FALSE)
  }
  .check_direction_orientation(sensor_rotation, azimuth_rotation, elevation)
  angles <- terra::values(g)
  .check_vector(angles[, 1], "numeric", allow_na = TRUE,
                name = "off_axis")
  .check_vector(angles[, 2], "numeric", allow_na = TRUE,
                name = "internal_azimuth")
  if (any(angles[, 1] < 0 | angles[, 1] > 180, na.rm = TRUE) ||
      any(angles[, 2] < 0 | angles[, 2] >= 360, na.rm = TRUE)) {
    stop("Internal angles must be in [0, 180] for `off_axis` and ",
         "[0, 360) for `internal_azimuth`.", call. = FALSE)
  }

  external_angles <- .external_direction_angles(
    angles * pi/180, sensor_rotation, azimuth_rotation, elevation) * 180/pi
  external_angles[!is.na(external_angles[, 2]) & external_angles[, 2] >= 360, 2] <- 0
  external <- terra::rast(g)
  terra::values(external) <- external_angles
  names(external) <- c("zenith", "azimuth")
  external
}

.check_direction_orientation <- function(sensor_rotation, azimuth_rotation,
                                          elevation) {
  .check_vector(sensor_rotation, "numeric", 1)
  .check_vector(azimuth_rotation, "numeric", 1)
  .check_vector(elevation, "numeric", 1)
  if (elevation < -90 || elevation > 90) {
    stop("`elevation` must lie within [-90, 90] degrees.", call. = FALSE)
  }
  invisible(TRUE)
}
