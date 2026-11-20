#' Interpolate in planar space
#'
#' @description
#' Interpolate values from canopy photographs using inverse distance weighting
#' (IDW) with k-nearest neighbors in image (planar) coordinates.
#' A radius limits neighbor search.
#'
#' @details
#' Defaults follow \insertCite{Lang2013;textual}{rcaiman}. Note that
#' `search_radius` is given in pixels but intended to approximate 15–20 deg in
#' angular terms. Therefore, this value needs fine-tuning based on image
#' resolution and radial projection. The proper workflow requires geometric
#' standardization with [reproject_image()] and [equidistant_coef]. For best
#' results, the interpolated quantity should be linearly related to scene
#' radiance--; for JPEG images, consider
#' [invert_gamma_correction()].
#'
#' @param sampling_coord `data.frame` with numeric columns `"x"` and
#'   `"y_down"` containing continuous image coordinates, and a numeric column
#'   with values to interpolate. Coordinates follow [y_down_from_y()]. Convert
#'   cell indices explicitly with [coord_from_indices()] when sampling cell
#'   centers. Subpixel positions are used without snapping to cells. Rows with
#'   missing or nonfinite coordinates or values are omitted. With no available
#'   samples, the result contains only `NA`.
#' @param col_id numeric or character vector of length one. The name or position
#'   of the column in `sampling_coord` containing the values to interpolate.
#' @param r numeric [terra::SpatRaster-class] with one layer. Image frame to
#'   which `sampling_coord` refers, with geometry set by [set_rcaiman_geometry()].
#'   Used only as a geometric template; cell values are ignored. Interpolation
#'   is evaluated at cell centers.
#' @param k numeric vector of length one. Must be positive and interpreted as
#'   integer. Number of nearest neighbors to retrieve.
#' @param p numeric vector of length one. Inverse distance weighting exponent.
#' @param search_radius numeric vector of length one. Maximum distance (pixels) to
#'   include neighbors.
#' @param rigorous logical vector of length one. By default (`FALSE`), the
#'   method interpolates as long as there is at least one point not farther than
#'   `search_radius`. If it is `TRUE`, it requires `k` points at range.
#' @param engine character vector of length one. Interpolation engine to be
#'   used. Allowed choices are `"cpp"` (default), an in-package C++ implementation,
#'   or `"terra"`, which delegates to [terra::interpIDW()], passing `k`, `p`,
#'   and `search_radius` to `maxPoints`, `power`, and `radius`, respectively,
#'   and when `rigorous = TRUE`, changes from `minPoints = 1` (default) to
#'   `minPoints = k`. When `"terra"` is used, arguments `parallel`, `cores`,
#'   `logical` and `leave_free` are ignored. Both engines use continuous input
#'   positions.
#'
#' @inheritParams apply_locally
#'
#' @note
#' Interpolation is performed in the image plane assuming Euclidean
#' geometry. When input points represent directions on the hemisphere, results
#' depend on the projection used to map angles to the plane. Thus, for
#' hemispherical data, the use of an equisolid-angle projection (see
#' [equal_area_coef]) is recommended.
#'
#' When interpolation is carried out using a dense set of points, errors
#' introduced by planar interpolation are expected to be negligible, i.e.,
#' increasing point density should not change the interpolated surface in a
#' meaningful way.
#'
#' No consistency checks are performed to ensure that `sampling_coord` and `r`
#' are geometrically compatible. Incorrect combinations may lead to invalid
#' outputs.
#'
#' @return Numeric [terra::SpatRaster-class] with one layer and the same geometry
#'   as `r`.
#'
#' @references \insertAllCited{}
#'
#' @seealso [panoramic()] for an alternative pipeline involving a temporary shift of image plane
#'
#' @export
#'
#' @examples
#' r <- set_rcaiman_geometry(terra::rast(nrows = 5, ncols = 7))
#' sampling_coord <- coord_from_indices(
#'   col = c(1, 7, 1, 7, 2, 5, 3, 6),
#'   row = c(1, 1, 5, 5, 3, 2, 4, 4)
#' )
#' sampling_coord$dn <- c(2, 7, 4, 11, 3, 9, 6, 10)
#' interpolated <- interpolate_planar(
#'   sampling_coord, r = r, search_radius = 20, parallel = FALSE
#' )
#' plot(interpolated)
#' points(sampling_coord$x, y_from_y_down(sampling_coord$y_down, nrow(r)))
interpolate_planar <- function(sampling_coord,
                               col_id = "dn",
                               r,
                               k = 3,
                               p = 2,
                               search_radius = 200,
                               rigorous = FALSE,
                               engine = "cpp",
                               parallel = TRUE,
                               cores = NULL,
                               logical = TRUE,
                               leave_free = 0) {

  col_id <- .resolve_required_column(col_id, sampling_coord)
  .check_columns(sampling_coord,
                 required_columns = c("x", "y_down", col_id))
  # Only the template grid is needed; do not read its cell values.
  if (!inherits(r, "SpatRaster") || terra::nlyr(r) != 1) {
    stop("`r` must be a single-layer SpatRaster.", call. = FALSE)
  }
  .check_vector(k, "integerish", 1, sign = "positive")
  .check_vector(p, "numeric", 1, sign = "positive")
  .check_vector(search_radius, "numeric", 1, sign = "positive")
  .check_vector(rigorous, "logical", 1)
  .assert_choice(engine, c("terra", "cpp"))
  for (column in unique(c("x", "y_down", col_id))) {
    .check_vector(sampling_coord[[column]], "numeric", allow_na = TRUE,
                  allow_inf = TRUE, name = column)
  }
  available <- is.finite(sampling_coord$x) &
    is.finite(sampling_coord$y_down) & is.finite(sampling_coord[[col_id]])
  sampling_coord <- sampling_coord[available, , drop = FALSE]

  if (engine == "cpp") {
    .check_vector(parallel, "logical", 1)
    .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
    .check_vector(logical, "logical", 1)
    .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

    if (parallel) {
      cores <- .cores(cores, logical, leave_free)
    } else {
      cores <- 1
    }
  }

  if (nrow(sampling_coord) == 0) {
    return(terra::setValues(r, NA_real_))
  }

  if (engine == "terra") {
    sampling_xy <- cbind(
      x = sampling_coord$x,
      y = y_from_y_down(sampling_coord$y_down, nrow(r))
    )
    v <- terra::vect(sampling_xy, crs = crs(r),
                      atts = data.frame(z = sampling_coord[, col_id]))

    ir <- terra::interpIDW(r, v, field = "z", radius = search_radius, power = p,
                           maxPoints = k, minPoints = if (rigorous) k else 1)
  }
  if (engine == "cpp") {

    target_indices <- terra::rowColFromCell(r, 1:terra::ncell(r))
    target_coord <- coord_from_indices(
      col = target_indices[, 2], row = target_indices[, 1]
    )

    res_vals <- cpp_idw_planar_kdtree(
      # Transition (2026-09-04): target_xy is the native engine's legacy
      # argument name. Both inputs here use coord, not the xy representation.
      # Keep its Euclidean-distance algorithm; rename that interface separately.
      target_xy = as.matrix(target_coord),
      obs_x = sampling_coord$x,
      obs_y = sampling_coord$y_down,
      obs_z = sampling_coord[[col_id]],
      k = k,
      p = p,
      rmax = search_radius,
      rigorous = rigorous,
      cores = cores
    )

    ir <- r
    ir[] <- res_vals
  }
  ir
}
