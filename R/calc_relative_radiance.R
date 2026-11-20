#' Normalize digital numbers by zenith value
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Compute relative radiance at sampled sky points using the zenith digital
#' number (DN) as reference.
#'
#' @details
#' Estimate zenith DN using inverse-distance weighting of the k nearest sky
#' points (weights proportional to \eqn{1 / d^2}, where \eqn{d} is angular
#' distance, see [apply_locally()]).
#' If selected points lie exactly at the zenith, only their DN values
#' contribute to the reference, with equal weights.
#'
#' @param sky_points `data.frame` with columns `"row"`, `"col"`, `"zenith"`, and
#'   `"azimuth"` (external angles in radians) as in the output of [add_directions()], and at least one
#'   additional column identified by argument `col_id`. Other columns are
#'   preserved.
#' @param col_id character or integer vector of length one. Name or position of
#'   the column in `sky_points` containing the DN values. DN is expected to be
#'   linearly related to radiance \insertCite{Lang2010}{rcaiman}.
#' @param k numeric vector of length one. Must be positive and interpreted as
#'   integer. Number of near-zenith points used to estimate the zenith DN using
#'   inverse distance weighting (power = 2). Ignored when `zenith_dn` is provided.
#' @param zenith_dn optional numeric vector of length one. DN at the zenith.
#'   Default `NULL`.
#'
#' @return List with named elements:
#' \describe{
#'   \item{`zenith_dn`}{A numeric vector of length one. Estimated DN at the zenith.}
#'   \item{`sky_points`}{A `data.frame` same as `sky_points` with additional
#'    column `relative_radiance` (relative radiance).
#'   }
#' }
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' sky_coef <- as.numeric(cie_table[4, 1:5])
#'
#' sun_angles <- data.frame(zenith = 45 * pi/180, azimuth = 10 * pi/180)
#'
#' m <- !is.na(d$zenith) & d$zenith < 80
#' sky <- cie_image(d, sun_angles, sky_coef, m = m)
#' plot(sky)
#' r <- sky * 20
#' names(r) <- "dn"
#'
#' directions <- fibonacci_directions(5)
#' masked_g <- g
#' masked_g[!m] <- NA
#' fov <- extract_fov(masked_g)
#' directions <- select_directions_in_fov(directions, fov)
#'
#' sky_points <- coord_from_internal_directions(directions, geometry) %>% add_indices()
#'
#' sky_points <- extract_dn(r, sky_points, use_window = FALSE)
#' sky_points <- add_directions(sky_points, d = d)
#' head(sky_points)
#' calc_relative_radiance(sky_points, "dn")
#'
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/DSCN4248.tif", package = "rcaiman")
#' metadata <- read.delim(gsub("tif", "tsv", path))
#' sensor_data <- read_caim(path)
#'
#' # Get specs
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' raw <- get_signature(my_entry, "raw")
#' cygm <- get_band_spec(raw, "cygm")
#' raw_geometry <- get_geometry_spec(cygm, "simple_method")
#' raw_radiometry <- get_radiometry_spec(raw_geometry, "photometric_sphere")
#'
#' # Represent image geometry
#' g <- set_internal_geometry(raw$dim, raw_geometry)
#' d <- direction_image(g)
#' m <- !is.na(g$off_axis)
#' w <- area_weight_image(g, raw_geometry)
#'
#' # Read and preprocess the image
#' caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
#' caim
#' r <- caim$Cyan
#' ff <- trend_surface_image(r, raw_radiometry$parameters[[metadata$aperture]]$Cyan, m = m)
#' r <- apply_flat_field(r, ff)
#'
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read_points(path) %>% add_indices()
#'
#' sky_points <- extract_dn(r, sky_points, use_window = TRUE)
#' sky_points <- add_directions(sky_points, d = d)
#' head(sky_points)
#' rr <- calc_relative_radiance(sky_points, "Cyan")
#' rr
#' }
calc_relative_radiance <- function(sky_points,
                                   col_id,
                                   k = 3,
                                   zenith_dn = NULL) {
  col_id <- .resolve_required_column(col_id, sky_points)
  .check_columns(
    sky_points,
    required_columns = c("row", "col", "zenith", "azimuth", col_id),
    reserved_columns = "relative_radiance"
  )
  .check_vector(zenith_dn, "numeric", 1, sign = "positive", allow_null = TRUE)

  if (is.null(zenith_dn)) {
    .check_vector(k, "integerish", 1, sign = "positive")
    if (nrow(sky_points) < k) {
      stop("`sky_points` must contain at least `k` rows.")
    }
    # Estimate zenith DN
    zenith_point <- data.frame(zenith = 0, azimuth = 0)
    nn <- nn_spherical(sky_points, zenith_point, k = k,
                       angle_columns = c("zenith", "azimuth"))

    zenith_dn <- apply_locally(
      sky_points,
      zenith_point,
      nn,
      fun = function(df, chi, ...) {
        b_k <- df[, col_id]
        w_k <- .inverse_distance_weights(chi)
        sum(w_k * b_k)
      },
      parallel = FALSE
    )$output

    if (min(sky_points$zenith) > .degree2radian(10))
      warning("The nearest sampled sky point lies at a zenith distance of ",
               round(.radian2degree(min(sky_points$zenith)), 2), " deg from the zenith." )
  }

  # Calculate relative radiance
  sky_points$relative_radiance <- sky_points[, col_id] / zenith_dn

  list(zenith_dn = zenith_dn,
       sky_points = sky_points)
}
