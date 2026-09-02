#' Triangulation-based interpolation in spherical space
#'
#' @description
#' Generate values at query locations using interpolation based on a
#' spherical Delaunay triangulation constructed from sampling directions.
#'
#' @details
#' The method operates on the unit sphere by converting zenith–azimuth
#' coordinates into 3D Cartesian vectors and constructing a triangulation
#' from those directions. The triangulation is obtained by computing the
#' convex hull of the augmented point set formed by the original directions
#' and their antipodal counterparts, and retaining only faces composed of
#' original points. This establishes a spherical Delaunay triangulation through the
#' convex-hull duality of points on the unit sphere.
#'
#' For each query direction, the enclosing spherical triangle is identified
#' using oriented edge tests. In boundary cases where multiple triangles are
#' valid due to numerical tolerance, the triangle whose normal vector is most
#' aligned with the query direction is selected, ensuring deterministic
#' behavior under rotation of the hemispherical domain.
#'
#' Interpolation is performed using barycentric weights derived from spherical
#' triangle areas (solid angles), yielding a continuous interpolated field over
#' the triangulated sphere.
#'
#' @inheritParams interpolate_spherical
#' @inheritParams interpolate_planar
#'
#' @return
#' `data.frame` equal to `query_points` with added columns `col_id` (by name) and
#' `"triangle_id"`, and without points outside the interpolation domian, i.e.,
#' for which `"triangle_id"` equals to `NA`.
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
#' query_indices <- fibonacci_indices(z, a, 1)
#' pts <- triangulate_spherical(sky_points, query_indices, "Cyan")
#'
#' sky <- interpolate_planar(pts, "Cyan", r,
#'                           k = 3, p = 2, search_radius = 6, rigorous = TRUE)
#' plot(sky)
#'
#' # Extrapolate and merge
#' # x11()
#' # plot(r)
#' # sun_angles <- click(c(z, a), 1) %>% as.numeric()
#' sun_angles <- c(zen_deg = 49.5, azi_deg = 27.4) #taken with above lines then hardcoded
#'
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
#' pts2 <- interpolate_spherical(pts, query_points, "Cyan", "dn_cie")
#'
#' sky <- interpolate_planar(pts2, "Cyan", r, k = 3, p = 2,
#'                           search_radius = 6, rigorous = TRUE)
#' plot(sky)
#' plot(r/sky)
#' }
triangulate_spherical <- function(sky_points,
                                  query_points,
                                  col_id) {

  col_id <- .resolve_required_column(col_id, sky_points)
  .check_columns(sky_points,
                 required_columns = c("row", "col", "zen", "azi", col_id)
  )

  .check_columns(query_points,
                 required_columns = c("row", "col", "zen", "azi"),
                 reserved_columns = c(col_id, "triangle_id")
  )

  .cross <- function(a, b) {
    c(
      a[2]*b[3] - a[3]*b[2],
      a[3]*b[1] - a[1]*b[3],
      a[1]*b[2] - a[2]*b[1]
    )
  }

  .dot <- function(a, b) sum(a * b)

  # área esférica (ángulo sólido)
  area <- function(a, b, c) {
    det_val <- abs(det(rbind(a, b, c)))
    dot_sum <- 1 + sum(a*b) + sum(b*c) + sum(c*a)
    2 * atan2(det_val, dot_sum)
  }

  # --- 1. datos ---
  pts3d <- as.matrix(.directions_to_cartesian(sky_points))
  vals  <- sky_points[[col_id]]
  n_orig <- nrow(pts3d)

  # --- 2. completar esfera ---
  pts3d_full <- rbind(pts3d, -pts3d)

  # --- 3. hull ---
  tri_full <- geometry::convhulln(pts3d_full, options = "Qt")

  # --- 4. filtrar triángulos válidos ---
  keep <- apply(tri_full, 1, function(t) all(t <= n_orig))
  tri <- tri_full[keep, , drop = FALSE]

  n_tri <- nrow(tri)

  # --- 5. precomputar aristas orientadas ---
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

  # --- 6. query points ---
  pix3d <- as.matrix(.directions_to_cartesian(query_points))

  tri_id <- rep(NA_integer_, nrow(pix3d))
  out <- rep(NA_real_, nrow(pix3d))
  eps <- 1e-12

  # 6.1. soluciona arbitrariedad en pertenencia para puntos en el borde
  tri_normals <- matrix(NA_real_, nrow = n_tri, ncol = 3)

  for (i in seq_len(n_tri)) {
    v1 <- pts3d[tri[i,1], ]
    v2 <- pts3d[tri[i,2], ]
    v3 <- pts3d[tri[i,3], ]

    n <- .cross(v2 - v1, v3 - v1)
    n <- n / sqrt(sum(n^2))  # normalizar

    tri_normals[i, ] <- n
  }

  # --- 7. loop ---
  for (i in seq_len(nrow(pix3d))) {
    p <- pix3d[i, ]

    s1 <- edges[,1,1]*p[1] + edges[,1,2]*p[2] + edges[,1,3]*p[3]
    s2 <- edges[,2,1]*p[1] + edges[,2,2]*p[2] + edges[,2,3]*p[3]
    s3 <- edges[,3,1]*p[1] + edges[,3,2]*p[2] + edges[,3,3]*p[3]

    idx <- which(s1 >= -eps & s2 >= -eps & s3 >= -eps)

    if (length(idx) == 0) next

    if (length(idx) == 1) {
      t_id <- idx
    } else {
      # elegir el triángulo cuya normal esté más alineada con p
      scores <- tri_normals[idx, , drop = FALSE] %*% p
      t_id <- idx[which.max(scores)]
    }

    t <- tri[t_id, ]
    tri_id[i] <- t_id

    v1 <- pts3d[t[1], ]
    v2 <- pts3d[t[2], ]
    v3 <- pts3d[t[3], ]

    A <- area(v1, v2, v3)
    if (A < 1e-14) next

    w1 <- area(p, v2, v3) / A
    w2 <- area(v1, p, v3) / A
    w3 <- area(v1, v2, p) / A

    out[i] <- w1 * vals[t[1]] +
      w2 * vals[t[2]] +
      w3 * vals[t[3]]
  }

  # --- salida ---
  query_points[[col_id]] <- out
  query_points$triangle_id <- tri_id
  query_points[!is.na(tri_id),]
}
