#' Interpolate in spherical space in point-wise manner
#'
#' @description
#' Generate new values on query locations using inverse distance weighting
#' (IDW) with k-nearest neighbors, computing distances in spherical coordinates.
#' Optionally blend with complementary source to fill voids.
#'
#' @details
#' The function is a wrapper of [apply_locally()]. When filling_source is
#' provided, interpolation follows the blending strategy proposed by
#' \insertCite{Lang2010;textual}{rcaiman}, where local inverse-distance weighted
#' estimates are combined with a complementary surface.
#'
#' \deqn{
#' B(\theta, \phi) =
#'   \frac{
#'     \sum_{k=1}^{n} w_k \, b_k + W \, M(\theta, \phi)
#'   }{
#'     \sum_{k=1}^{n} w_k + W
#'   }
#' }
#'
#' with \eqn{w_k = \frac{1}{\chi_k^p}}
#'
#' where \eqn{B(\theta, \phi)} is the value to interpolate in direction \eqn{(\theta,
#' \phi)}, calculated with \eqn{k = 0, 1,... n} closest values \eqn{b_k}, and
#' the filling source \eqn{M(\theta, \phi)}, with user defined weigth \eqn{W}
#' here equal to `w`.
#'
#' The original version propose a normalization of \eqn{w_k} based on the
#' available neighborhood. However, it may produce context-dependent behavior. For
#' example, with `k = 1`, local normalization would force the single neighbor to
#' full weight regardless of its actual angular distance, making the blend
#' insensitive to proximity.
#'
#' To preserve a stable interpretation of `w` across interpolation contexts, this
#' implementation normalizes inverse-distance weights using fixed angular bounds
#' defined by the interpolation domain rather than by the observed neighborhood
#' itself. Specifically, weights are computed from angular distances in
#' spherical space and min-max rescaling using:
#'
#' \deqn{
#'  w_{min} = 1 / chi_{max}^p
#'  }
#' \deqn{
#'  w_{max} = 1 / chi_{min}^p
#'  }
#'
#' where, \eqn{chi_{max}} is equal to `search_radius`  and \eqn{chi_{min}} is
#' equal to `min_dist`, both internally transformed to radians. Rescaling is
#' carried out with [normalize_min_max()] with `clip = TRUE`.
#'
#' This yields bounded weights with stable meaning across query locations while
#' preserving the expected transition between local dominance and model-based
#' filling.
#'
#' For efficiency, interpolation runs in point-wise manner, meaning that
#' computing workload can regulated with the number of query points. To obtain a
#' raster use interpolate_planar() with the query points as input. Use
#' `raster_grid_in_circle(ncol(r), 1)` to produce one query point per pixel.
#'
#' @param sky_points `data.frame` with columns `"zen"` and `"azi"` (see
#'   [add_directions()]), and one additional numeric column with values to
#'   interpolate, identified by argument `col_id`. Additional columns are ignored.
#' @param query_points `data.frame` defining locations where interpolation is
#'   evaluated. Must contain columns `""zen""` and `""azi""` (see
#'   [add_directions()]). Cannot contain a column named as `col_id` (by name).
#'   Additional columns are ignored.
#' @param filling_source numeric or character vector of length one. The name or
#'   position of the column in `query_points` providing complementary values.
#' @param min_dist numeric vector of length one. Minimum angular distance in
#'   degrees considered for min-max rescaling. See *Details*.
#' @param w numeric vector of length one or single-layer
#'   [terra::SpatRaster-class]. Weight assigned to `filling_source` in the blend
#'   with local estimates (see Eq. 6 in \insertCite{Lang2010;textual}{rcaiman}).
#' @param rule character vector of length one. Either `"any"` or `"all"`.
#'   With `"any"`, pixels within `search_radius` of at least one sample are
#'   interpolated. With `"all"`, pixels are interpolated only if the `k`
#'   nearest neighbors lie within `search_radius`. If `k = 1`, both are
#'   equivalent.
#'
#' @inheritParams apply_locally
#' @inheritParams interpolate_planar
#'
#' @references \insertAllCited{}
#'
#' @return
#' `data.frame` equal to `query_points` with added columns `col_id` (by name)
#' and without points outside the interpolation domian, i.e., for which
#' `col_id` equals to `NA`.
#'
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
#' rspec <- get_radiometry_spec(spec, "simple_method")
#'
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#'
#' # See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#'
#' query_points <- raster_grid_in_circle(z, a, 10)
#'
#' pts <- interpolate_spherical(
#'   sky_points,
#'   query_points,
#'   col_id = "Cyan",
#'   k = 3,
#'   p = 2,
#'   search_radius = 20,
#'   rule = "any"
#' )
#'
#' head(pts)
#' sky <- interpolate_planar(pts, "Cyan", r,
#'                           k = 3, p = 2, search_radius = 13, rigorous = TRUE)
#' plot(sky)
#'
#' # x11()
#' # plot(r)
#' # sun_angles <- click(c(z, a), 1) %>% as.numeric()
#' sun_angles <- c(zen_deg = 49.5, azi_deg = 27.4) #taken with above lines then hardcoded
#' rr <- calc_relative_radiance(sky_points, 3)
#'
#' set.seed(7)
#' model <- fit_cie_model(rr, sun_angles,
#'                        general_sky_type = "Clear")
#'
#' sky_cie <- cie_image(z, a, model$sun_angles, model$coef)
#' sky_cie <- sky_cie * rr$zenith_dn
#' names(sky_cie) <- "dn_cie"
#'
#' query_points <- raster_grid_in_circle(z, a, 5)
#' query_points <- extract_dn(sky_cie, query_points, use_window = FALSE)
#' head(query_points)
#'
#' pts <- interpolate_spherical(
#'   rr$sky_points,
#'   query_points,
#'   col_id = "Cyan",
#'   filling_source = "dn_cie",
#'   min_dist = 1,
#'   w = 0.01,
#'   k = 3,
#'   p = 2,
#'   search_radius = 20,
#'   rule = "any"
#' )
#'
#' head(pts)
#' sky <- interpolate_planar(pts, "Cyan", r,
#'                           k = 3, p = 2, search_radius = 7, rigorous = TRUE)
#' plot(sky)
#' plot(r/sky)
#' }
interpolate_spherical <- function(sky_points,
                                  query_points,
                                  col_id,
                                  filling_source = NULL,
                                  w = 1,
                                  k = 3,
                                  p = 2,
                                  search_radius = 20,
                                  min_dist = 1,
                                  rule = "any",
                                  parallel = FALSE,
                                  cores = NULL,
                                  logical = TRUE,
                                  leave_free = 1
                              ) {

  col_id <- .resolve_required_column(col_id, sky_points)
  .check_columns(sky_points,
                 required_columns = c("row", "col", "zen", "azi", col_id),
                 reserved_columns = NULL
                 )

  if (is.null(filling_source)) {
    .check_columns(query_points,
                   required_columns = c("row", "col", "zen", "azi"),
                   reserved_columns = col_id
                   )
  } else {
    filling_source <- .resolve_required_column(filling_source, query_points)
    .check_columns(query_points,
                   required_columns = c("row", "col", "zen", "azi", filling_source),
                   reserved_columns = col_id
                   )
    .check_vector(min_dist, "numeric", 1, sign = "positive")
  }

  .check_vector(w, "numeric", 1, sign = "nonnegative")
  .check_vector(p, "numeric", 1, sign = "positive")

  neighbors <- nn_spherical(sky_points,
                            query_points,
                            k = k,
                            search_radius = search_radius,
                            rule = "any")

  if (is.null(filling_source)) {
    eval_points <- apply_locally(
      sky_points,
      query_points,
      neighbors,
      fun = function(df, chi, ...) {
        b_k <- df[, col_id]
        w_k <- 1 / chi^p
        sum(w_k * b_k) / sum(w_k)
      },
      parallel = parallel,
      cores = cores,
      logical = logical,
      leave_free = leave_free
    )
  } else {
    # eq6 from Lang et al. 2010
    .eq6 <- function(w_k, b_k, W, M) {
      (sum(w_k * b_k) + W * M) / (sum(w_k) + W)
    }
    chi_mx <- search_radius * pi/180
    chi_mn <- min_dist * pi/180

    eval_points <- apply_locally(
      sky_points,
      query_points,
      neighbors,
      fun = function(df, chi, i, ...) {
        M <- query_points[i, filling_source]
        if (length(chi) != 0) {
          w_k <- 1 / chi^p
          w_mn <- 1 / chi_mx^p
          w_mx <- 1 / chi_mn^p
          w_k <- normalize_min_max(w_k, w_mn, w_mx, clip = TRUE)
          .eq6(w_k, df[, col_id], w, M)
        } else {
          M
        }
      },
      parallel = parallel,
      cores = cores,
      logical = logical,
      leave_free = leave_free
    )
  }

  i <- !is.na(eval_points$output)
  eval_points <- eval_points[i,]
  names(eval_points)[ncol(eval_points)] <- col_id
  eval_points
}
