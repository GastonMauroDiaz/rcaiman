#' Normalize digital numbers by zenith value
#'
#' @description
#' Compute relative radiance at sampled sky points using the zenith digital
#' number (DN) as reference.
#'
#' @details
#' Estimate zenith DN using inverse-distance weighting of the k nearest sky
#' points (weights proportional to \eqn{1 / d^2}, where \eqn{d} is angular
#' distance, see [apply_locally()]).
#'
#' @param sky_points `data.frame` with columns `"row"`, `"col"`, `"zen"`, and
#'   `"azi"` as in the output of [add_directions()], and at least one
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
#'    column `rr` (relative radiance).
#'   }
#' }
#'
#' @references \insertAllCited{}
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
#'
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' display_caim(r, sampling_points = sky_points)
#' sky_points <- extract_dn(r, sky_points, use_window = TRUE)
#' sky_points <- add_directions(sky_points, z, a)
#' calc_relative_radiance(sky_points, "Cyan")
#' }
calc_relative_radiance <- function(sky_points,
                                   col_id,
                                   k = 3,
                                   zenith_dn = NULL) {
  col_id <- .resolve_required_column(col_id, sky_points)
  .check_columns(
    sky_points,
    required_columns = c("row", "col", "zen", "azi", col_id),
    reserved_columns = "rr"
  )
  .check_vector(zenith_dn, "numeric", 1, sign = "positive", allow_null = TRUE)

  if (is.null(zenith_dn)) {
    .check_vector(k, "integerish", 1, sign = "positive")
    if (nrow(sky_points) < k) {
      stop("`sky_points` must contain at least `k` rows.")
    }
    # Estimate zenith DN
    zenith_point <- data.frame(zen = 0, azi = 0)
    nn <- nn_spherical(sky_points, zenith_point, k = k)

    zenith_dn <- apply_locally(
      sky_points,
      zenith_point,
      nn,
      fun = function(df, chi, ...) {
        b_k <- df[, col_id]
        w_k <- 1 / chi^2
        sum(w_k * b_k) / sum(w_k)
      },
      parallel = FALSE
    )$output

    if (min(sky_points$zen) > .degree2radian(10))
      warning("The nearest sampled sky point lies at a zenith distance of ",
               round(.radian2degree(min(sky_points$zen)), 2), " deg from the zenith." )
  }

  # Calculate relative radiance
  sky_points$rr <- sky_points[, col_id] / zenith_dn

  list(zenith_dn = zenith_dn,
       sky_points = sky_points)
}
