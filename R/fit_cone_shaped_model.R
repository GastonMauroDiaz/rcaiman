#' Fit cone-shaped model
#'
#' @description
#' Fit a polynomial model to predict sky digital numbers from directions (see
#' [add_directions()]) using data sampled from a canopy photograph.
#'
#' @details
#' Depending on `method`, it can fit:
#' \describe{
#'   \item{*A zenith-only quadratic model*}{\deqn{sDN = a + b\theta + c\theta^2}}
#'   \item{*A zenith-plus-azimuth model, adding sinusoidal terms*}{\deqn{sDN = a + b\theta + c\theta^2 + d\sin(\phi) + e\cos(\phi)}}
#' }
#'
#' See \insertCite{Diaz2018;textual}{rcaiman} for additional details.
#'
#' See [sky_via_local_models()] for an extension of this approach that fits
#' multiple local cone-shaped models across hemispherical space.  By combining
#' local predictions, that approach can capture complex sky radiance patterns
#' without restricting the reconstructed sky surface to a single cone-shaped
#' structure.
#'
#' @param col_id character or integer vector of length one. Name or position of
#'   the column in `sky_points` containing the values sampled from the image.
#' @param method character. Model type to fit:
#'   \describe{
#'     \item{`"zenith_only"`}{Quadratic polynomial in zenith angle.}
#'     \item{`"zenith_n_azimuth"`}{Quadratic polynomial in zenith plus
#'       sinusoidal terms in azimuth.}
#'   }
#'
#' @inheritParams calc_relative_radiance
#'
#' @return List with the following components:
#' \describe{
#'   \item{`fun`}{Function taking `zen` and `azi` (angle units can be declared
#'     via the argument `units`) and returning predicted values.}
#'   \item{`model`}{`lm` object fitted by [stats::lm()].}
#' }
#' Returns `NULL` (with a warning) if the number of input points is fewer than 20.
#'
#' @export
#'
#' @references \insertAllCited{}
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
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' plot(r)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#'
#' model <- fit_cone_shaped_model(sky_points, col_id = "Cyan")
#' summary(model$model)
#' sky_cs <- model$fun(z, a)
#' sky_cs
#' plot(sky_cs)
#'
#' z_mini <- zenith_image(50, equidistant_coef)
#' sky_cs <- model$fun(z_mini, azimuth_image(z_mini))
#' persp(sky_cs, theta = 90, phi = 20)
#'
#' }
fit_cone_shaped_model <- function(sky_points,
                                  col_id,
                                  method = "zenith_n_azimuth") {
  col_id <- .resolve_required_column(col_id, sky_points)
  .check_columns(sky_points,
                 required_columns = c("col", "row", "zen", "azi", col_id),
                 reserved_columns = NULL)
  .assert_choice(method, c("zenith_only", "zenith_n_azimuth"))

  Blue <- sky_points[, col_id]
  Zenith <- sky_points$zen
  Azimuth <- sky_points$azi

  if (length(Blue) >= 20) {
    if (method == "zenith_n_azimuth") {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE) +
        sin(Azimuth) + cos(Azimuth))

      # Only to avoid note from check, code is OK without this line.
      a <- b <- d <- e <- NA

      .sky_fun <- function(zen, azi, units = "degree") {
        .assert_choice(units, c("degree", "radian"))
        if (units == "degree") {
          zen <- .degree2radian(zen)
          azi <- .degree2radian(azi)
        }
        x <- coefficients(model)
        x[is.na(x)] <- 0
        for (i in 1:5) assign(letters[i], x[i])
        a + b * zen + c * zen^2 +
          d * sin(azi) + e * cos(azi)
      }
    } else {
      model <- lm(Blue ~ poly(Zenith, 2, raw = TRUE))
      .sky_fun <- function(zenith, azimuth) {
        zenith <- .degree2radian(zenith)
        azimuth <- .degree2radian(azimuth)
        x <- coefficients(model)
        x[is.na(x)] <- 0
        for (i in 1:5) assign(letters[i], x[i])
        a + b * zenith + c * zenith^2
      }
    }
    return(list(fun = .sky_fun, model = model))
  } else {
    warning("Insufficient number of points to attempt model fitting")
    return(NULL)
  }
}
