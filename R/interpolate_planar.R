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
#' standarization with [reproject_image()] and [equidistant_coef]. For best
#' results, the interpolated quantity should be linearly related to scene
#' radiance--; for JPEG images, consider
#' [invert_gamma_correction()].
#'
#' @param sampling_points `data.frame` with columns `"row"`, `"col"`, and one additional
#'   numeric column with values to interpolate. Typically returned by
#'   [extract_dn()].
#' @param col_id numeric or character vector of length one. The name or position
#'   of the column in `sampling_points` containing the values to interpolate.
#' @param r numeric [terra::SpatRaster-class] with one layer. Image from which
#'   `sampling_points` were derived, or another raster with the same number of rows
#'   and columns. Used only as geometric template; cell values are ignored.
#' @param k numeric vector of length one. Must be positive and interpreted as
#'   integer. Number of nearest neighbors to retrieve.
#' @param p numeric vector of length one. Inverse distance weighting exponent.
#' @param search_radius numeric vector of length one. Maximum distance (pixels) to
#'   include neighbors.
#' @param rigorous logical vector of length one. By default (`FALSE`), the
#'   method interpolates as long as there is at least one point not farther than
#'   `search_radius`. If it is `TRUE`, it requires `k` points at range.
#' @param engine character vector of length one. Interpolation engine to be
#'   used. Allowed choises are `"cpp"` (default), an in-package C++ implementation,
#'   or `"terra"`, which delegates to [terra::interpIDW()], passing `k`, `p`,
#'   and `search_radius` to `maxPoints`, `power`, and `radius`, respectively,
#'   and when `rigorous = TRUE`, changes from `minPoints = 1` (default) to
#'   `minPoints = k`. When `"terra"` is used, arguments `parallel`, `cores`,
#'   `logical` and `leave_free` are ignored. The effective search neighborhood
#'   differs slightly from engines, likely due to continuous vs rasterized
#'   distance evaluation.
#'
#' @inheritParams apply_locally
#'
#' @note
#' Interpolation is performed in the image plane assuming Euclidean
#' geometry. When input points represent directions on the hemisphere,results
#' depend on the projection used to map angles to the plane. Thus, for
#' hemispherical data, the use of an equisolid-angle projection (see
#' [equal_area_coef]) is recommended.
#'
#' When interpolation is carried out using a dense set of points, errors
#' introduced by planar interpolation are expected to be negligible, i.e.,
#' increasing point density should not change the interpolated surface in a
#' meaninful way.
#'
#' No consistency checks are performed to ensure that `sampling_points` and `r`
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
#'
#' sky <- interpolate_planar(sky_points, r, col_id = 3)
#' plot(sky)
#' plot(r/sky)
#' }
interpolate_planar <- function(sampling_points,
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

  col_id <- .resolve_required_column(col_id, sampling_points)
  .check_columns(sampling_points, required_columns = c("row", "col", col_id))
  .assert_single_layer(r)
  .check_vector(k, "integerish", 1, sign = "positive")
  .check_vector(p, "numeric", 1, sign = "positive")
  .check_vector(search_radius, "numeric", 1, sign = "positive")
  .check_vector(rigorous, "logical", 1)
  .assert_choice(engine, c("terra", "cpp"))

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

  cells <- terra::cellFromRowCol(r, sampling_points$row, sampling_points$col)
  xy <- terra::xyFromCell(r, cells)

  if (engine == "terra") {
    v <- terra::vect(xy, crs = crs(r), atts = data.frame(z = sampling_points[, col_id]) )

    ir <- terra::interpIDW(r, v, field = "z", radius = search_radius, power = p,
                           maxPoints = k, minPoints = if (rigorous) k else 1)
  }
  if (engine == "cpp") {

    target_rc <- terra::rowColFromCell(r, 1:terra::ncell(r))

    res_vals <- cpp_idw_planar_kdtree(
      target_xy = as.matrix(target_rc[, 2:1]),
      obs_x = sampling_points$col,
      obs_y = sampling_points$row,
      obs_z = sampling_points[[col_id]],
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
