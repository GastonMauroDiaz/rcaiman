#' Triangulation-based interpolation in spherical space
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Interpolate sampled values at query directions using a triangulation on
#' the unit sphere and weights based on spherical triangle areas.
#'
#' @details
#' Both tables must use the same reference frame. Select external or internal
#' directions with `angle_columns`; angles are in radians, following
#' [add_directions()]. Cell indices and image coordinates are not needed.
#'
#' Directions are converted to 3D unit vectors. [geometry::convhulln()] with
#' option `"Qt"` constructs the convex hull of the observations and their
#' antipodal counterparts. Only triangular faces composed entirely of original
#' observations are retained. Their spherical triangles define the interpolation
#' domain. The augmented point set must support a three-dimensional convex hull;
#' errors from a degenerate configuration are passed through without perturbing
#' the directions.
#'
#' Oriented edge tests identify candidate triangles for each query, with a
#' numerical tolerance of \eqn{10^{-12}}. When several candidates are found,
#' the function chooses the triangle whose computed unit normal has the largest
#' dot product with the query direction. An exact tie selects the first candidate
#' in the returned triangulation order. Repeating the calculation with the same
#' input tables, row order, and software configuration gives the same result.
#' Triangle numbering is local to that triangulation, not a persistent spatial
#' identifier.
#'
#' For a query \eqn{p} in a triangle with vertices \eqn{a, b, c}, weights are
#' \eqn{A(p,b,c)/A(a,b,c)}, \eqn{A(a,p,c)/A(a,b,c)}, and
#' \eqn{A(a,b,p)/A(a,b,c)}, where \eqn{A} is spherical area in steradians.
#' The prediction is the weighted sum of the three observed values. Triangles
#' with area below \eqn{10^{-14}} steradians do not supply an estimate. Values
#' are not normalized or clipped. These area ratios do not in general give
#' identical edge values from adjacent triangles; continuity across their
#' boundaries is not guaranteed.
#'
#' Queries outside the triangulation are retained with `NA`. No extrapolation
#' or background filling is performed. For a separate interpolation or blending
#' step, see [interpolate_spherical()]. Use [add_coord()] when retained cell
#' indices need to supply [interpolate_planar()].
#'
#' @param sampling_directions data frame with numeric angular columns selected
#'   by `angle_columns` and a numeric value column selected by `col_id`. Samples
#'   with missing values or directions are excluded before triangulation.
#'   Infinite values are not allowed. Additional columns are ignored.
#' @param query_directions data frame with the same angular columns. All rows
#'   and additional columns are preserved. Must not already contain the
#'   sampled-value column or `triangle_id`.
#' @param col_id character or integer vector of length one. Name or position
#'   of the sampled-value column in `sampling_directions`. Its name must differ
#'   from the angular columns and `triangle_id`.
#' @param angle_columns character vector of length two, naming the polar and
#'   azimuthal columns in that order. Defaults to `c("zenith", "azimuth")`;
#'   for internal directions, use `c("off_axis", "internal_azimuth")`.
#'   Polar angles must lie in \eqn{[0, \pi]} and azimuths in
#'   \eqn{[0, 2\pi)}. A missing azimuth at polar angle zero or \eqn{\pi}
#'   still defines a valid direction. Other missing directions cannot be used.
#'
#' @return Data frame with the same rows, order, and columns as
#'   `query_directions`, plus the sampled-value column and integer
#'   `triangle_id`. Both added columns are `NA` for queries without an eligible
#'   triangle or with missing directions. With fewer than three usable samples,
#'   or no retained triangles, all estimates are `NA`.
#'
#' @seealso [add_directions()], [interpolate_spherical()], [interpolate_planar()]
#' @export
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 60)
#' m <- !is.na(d$zenith) & !is.na(d$azimuth) & d$zenith <= 90
#'
#' # Illustrative values; with measurements, select suitable sky observations.
#' r <- 1000 + 200 * (d$zenith * pi / 180)^2 +
#'   80 * sin(d$azimuth * pi / 180)
#' names(r) <- "Blue"
#' sampling_indices <- expand.grid(col = seq(1, 80, 5), row = seq(1, 60, 5))
#' sampling_indices <- rem_unselected(sampling_indices, m)
#' samples <- extract_dn(r, sampling_indices, use_window = FALSE)
#' samples <- add_directions(samples, d = d)
#'
#' query_indices <- expand.grid(col = seq(1, 80, 3), row = seq(1, 60, 3))
#' query_indices <- rem_unselected(query_indices, m)
#' queries <- add_directions(query_indices, d = d)
#' interpolated <- triangulate_spherical(samples, queries, "Blue")
#' head(interpolated)
#'
#' # The planar interpolation is a separate step on the image grid.
#' sky <- interpolate_planar(add_coord(interpolated), "Blue", r,
#'                           search_radius = 6, parallel = FALSE)
#' sky[!m] <- NA
#' plot(sky)
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
#' # Read cie model
#' path <- system.file("external/cie_model.RDS", package = "rcaiman")
#' model <- readRDS(path)
#'
#' sky_cie <- cie_image(d, model$sun_angles, model$coef) * model$relative_radiance$zenith_dn
#' names(sky_cie) <- "above-canopy dn"
#' plot(sky_cie)
#'
#' # Binarize
#' bin <- binarize_with_thr(r/sky_cie, 0.5)
#' bin <- distance_to_black(bin) > 2
#' plot(bin)
#'
#' # Dense sampling of the sky DN
#' seg <- healpix_segmentation(d, 32)
#' seg[!bin] <- 0
#' features <- compute_features(r, seg, d, w)
#' sky_points <- features[features$effective_pixel_area > 3, ]
#' sky_points <- internal_directions_from_external(sky_points) %>%
#'   coord_from_internal_directions(., raw_geometry) %>% add_indices()
#' head(sky_points)
#' nrow(sky_points)
#'
#' # Prepare interpolation
#' fov <- extract_fov(g)
#' query_directions <- fibonacci_directions(1)
#' query_directions <- select_directions_in_fov(query_directions, fov)
#' query_directions <- external_directions_from_internal(query_directions)
#'
#' # Interpolate without blending
#' interpolated <- triangulate_spherical(
#'   sky_points,
#'   query_directions,
#'   col_id = "max"
#' )
#' head(interpolated)
#' interpolated <- coord_from_internal_directions(interpolated, raw_geometry)
#' plot(interpolate_planar(interpolated, "max", r,  search_radius = 10))
#'
#' # Get complementary values for blending
#' query_directions <- coord_from_internal_directions(query_directions, raw_geometry) %>%
#'   add_indices()
#' query_directions <- extract_dn(sky_cie, query_directions, use_window = FALSE)
#'
#' # Interpolate and blend
#' blended <- interpolate_spherical(
#'   interpolated,
#'   query_directions,
#'   col_id = "max",
#'   filling_source = "above-canopy dn",
#'   filling_weight = 0.01,
#'   min_dist = 1,
#'   rule = "all",
#'   search_radius = 20
#' )
#' head(blended)
#'
#' sky <- interpolate_planar(blended, "max", r,  search_radius = 10)
#' plot(sky)
#' plot(r/sky)
#' }
triangulate_spherical <- function(sampling_directions,
                                  query_directions,
                                  col_id,
                                  angle_columns = c("zenith", "azimuth")) {
  .check_vector(angle_columns, "character", 2)
  if (anyDuplicated(angle_columns)) {
    stop("`angle_columns` must name two different columns.", call. = FALSE)
  }
  .check_columns(sampling_directions, required_columns = angle_columns)
  if (is.numeric(col_id)) .check_vector(col_id, "integerish", 1, sign = "positive")
  col_id <- .resolve_required_column(col_id, sampling_directions)
  if (col_id %in% c(angle_columns, "triangle_id")) {
    stop("The sampled-value column must differ from the angular columns and `triangle_id`.",
         call. = FALSE)
  }
  .check_columns(query_directions, required_columns = angle_columns,
                 reserved_columns = c(col_id, "triangle_id"))
  vals <- sampling_directions[[col_id]]
  .check_vector(vals, "numeric", allow_na = TRUE, name = col_id)
  pts3d <- as.matrix(.spherical_neighbor_vectors(sampling_directions, angle_columns))
  pix3d <- as.matrix(.spherical_neighbor_vectors(query_directions, angle_columns))
  available <- stats::complete.cases(pts3d) & !is.na(vals)
  pts3d <- pts3d[available, , drop = FALSE]
  vals <- vals[available]
  result <- query_directions
  result[[col_id]] <- rep(NA_real_, nrow(result))
  result$triangle_id <- rep(NA_integer_, nrow(result))
  if (nrow(pts3d) < 3L || !any(stats::complete.cases(pix3d))) return(result)

  .cross <- function(a, b) {
    c(
      a[2]*b[3] - a[3]*b[2],
      a[3]*b[1] - a[1]*b[3],
      a[1]*b[2] - a[2]*b[1]
    )
  }

  .dot <- function(a, b) sum(a * b)

  # Spherical triangle area (solid angle).
  area <- function(a, b, c) {
    det_val <- abs(det(rbind(a, b, c)))
    dot_sum <- 1 + sum(a*b) + sum(b*c) + sum(c*a)
    2 * atan2(det_val, dot_sum)
  }

  n_orig <- nrow(pts3d)

  # Augment with antipodal directions, retaining only original-point faces.
  pts3d_full <- rbind(pts3d, -pts3d)

  tri_full <- geometry::convhulln(pts3d_full, options = "Qt")

  keep <- apply(tri_full, 1, function(t) all(t <= n_orig))
  tri <- tri_full[keep, , drop = FALSE]

  n_tri <- nrow(tri)
  if (n_tri == 0L) return(result)

  # Precompute oriented edges.
  edges <- array(NA_real_, dim = c(n_tri, 3, 3))

  for (i in seq_len(n_tri)) {
    v1 <- pts3d[tri[i,1], ]
    v2 <- pts3d[tri[i,2], ]
    v3 <- pts3d[tri[i,3], ]

    e1 <- .cross(v1, v2)
    e2 <- .cross(v2, v3)
    e3 <- .cross(v3, v1)

    if (.dot(e1, v3) < 0) e1 <- -e1
    if (.dot(e2, v1) < 0) e2 <- -e2
    if (.dot(e3, v2) < 0) e3 <- -e3

    edges[i,1,] <- e1
    edges[i,2,] <- e2
    edges[i,3,] <- e3
  }

  tri_id <- rep(NA_integer_, nrow(pix3d))
  out <- rep(NA_real_, nrow(pix3d))
  eps <- 1e-12

  # Normals used to resolve multiple candidate triangles.
  tri_normals <- matrix(NA_real_, nrow = n_tri, ncol = 3)

  for (i in seq_len(n_tri)) {
    v1 <- pts3d[tri[i,1], ]
    v2 <- pts3d[tri[i,2], ]
    v3 <- pts3d[tri[i,3], ]

    n <- .cross(v2 - v1, v3 - v1)
    n <- n / sqrt(sum(n^2))

    tri_normals[i, ] <- n
  }

  for (i in which(stats::complete.cases(pix3d))) {
    p <- pix3d[i, ]

    s1 <- edges[,1,1]*p[1] + edges[,1,2]*p[2] + edges[,1,3]*p[3]
    s2 <- edges[,2,1]*p[1] + edges[,2,2]*p[2] + edges[,2,3]*p[3]
    s3 <- edges[,3,1]*p[1] + edges[,3,2]*p[2] + edges[,3,3]*p[3]

    idx <- which(s1 >= -eps & s2 >= -eps & s3 >= -eps)

    if (length(idx) == 0) next

    if (length(idx) == 1) {
      t_id <- idx
    } else {
      # Choose the most aligned normal; which.max breaks exact ties by order.
      scores <- tri_normals[idx, , drop = FALSE] %*% p
      t_id <- idx[which.max(scores)]
    }

    t <- tri[t_id, ]

    v1 <- pts3d[t[1], ]
    v2 <- pts3d[t[2], ]
    v3 <- pts3d[t[3], ]

    A <- area(v1, v2, v3)
    if (A < 1e-14) next
    tri_id[i] <- t_id

    w1 <- area(p, v2, v3) / A
    w2 <- area(v1, p, v3) / A
    w3 <- area(v1, v2, p) / A

    out[i] <- w1 * vals[t[1]] +
      w2 * vals[t[2]] +
      w3 * vals[t[3]]
  }

  result[[col_id]] <- out
  result$triangle_id <- tri_id
  result
}
