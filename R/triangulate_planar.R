#' Interpolate with Delaunay triangulation
#'
#' @description
#' Interpolate pixel values from irregularly spaced points using planar
#' surfaces defined by Delaunay triangles.
#'
#' @details
#' Each triangle defines a plane estimated from its three vertices, and
#' interpolation is performed by evaluating the plane equation at each
#' pixel position. This approach produces a continuous surface with
#' sharp transitions at triangle borders.
#'
#' The function internally uses [terra::delaunay()] to construct the
#' triangulated mesh.
#'
#' @note
#' Interpolation is performed in the image plane assuming Euclidean geometry.
#' When input points represent directions on the hemisphere, triangle edges do
#' not represent geodesics and therefore results depend on the projection used
#' to map angles to the plane. Thus, for hemispherical data, the use of an
#' equisolid-angle projection (see [equal_area_coef]) is recommended.
#'
#' When interpolation is carried out using a dense set of points, planar
#' triangulation can be interpreted as a local approximation to spherical
#' geometry. In this case, the angular spacing of query points should be small
#' relative to the spatial variability of the interpolated field. Under this
#' condition, discretization errors introduced by planar interpolation are
#' expected to be negligible, i.e., increasing point density should not change
#' the interpolated surface in a meaningful way.
#'
#' No consistency checks are performed to ensure that `sampling_coord` and `r`
#' are geometrically compatible. Incorrect combinations may lead to invalid
#' outputs.
#'
#' Vertex values are currently retrieved through an auxiliary raster. Distinct
#' positions within the same cell can therefore share a value. This existing
#' limitation is separate from the continuous-coordinate input contract.
#'
#' @param return character vector of length one. Either `"raster"`
#'   (default) to return only the interpolation result, or `"all"`
#'   to also return the Delaunay triangles and the convex hull.
#'
#' @inheritParams interpolate_planar
#'
#' @seealso [panoramic()] for an alternative pipeline involving a temporary shift of image plane
#'
#' @return
#' If `return = "raster"`, a numeric [terra::SpatRaster-class] with one layer
#' and the same geometry as `r`.
#'
#' If `return = "all"`, a list with named elements:
#' \describe{
#'   \item{`raster`}{the [terra::SpatRaster-class] described above.}
#'   \item{`mesh`}{a [terra::SpatVector-class] with the Delaunay triangles.}
#'   \item{`convex_hull`}{a [terra::SpatVector-class] with the convex hull.}
#' }
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
#' interpolated <- triangulate_planar(sampling_coord, r = r)
#' plot(interpolated)
#' points(sampling_coord$x, y_from_y_down(sampling_coord$y_down, nrow(r)))
triangulate_planar <- function(sampling_coord, col_id = "dn", r,
                                return = "raster") {
  col_id <- .resolve_required_column(col_id, sampling_coord)
  .check_columns(sampling_coord,
                 required_columns = c("x", "y_down", col_id))
  # Correction (2026-09-04): template cell values are ignored, including NA.
  .assert_single_layer(r, allow_na = TRUE)
  .assert_choice(return, c("raster", "all"))

  # Compute plane from three points using cross product
  # Points: numeric vectors of length 3 (x,y,z)
  .plane_from_points <- function(p1, p2, p3) {
    u <- p2 - p1
    v <- p3 - p1

    n <- c(
      u[2] * v[3] - u[3] * v[2],
      u[3] * v[1] - u[1] * v[3],
      u[1] * v[2] - u[2] * v[1]
    )

    n
    a <- n[1]
    b <- n[2]
    c <- n[3]
    d <- - sum(n * p1)  # d = -n . P1

    # z = alpha*x + beta*y + gamma
    unname(c(-a / c, -b / c, -d / c))
  }

  sampling_xy <- cbind(
    x = sampling_coord$x,
    y = y_from_y_down(sampling_coord$y_down, nrow(r))
  )

  p <- terra::vect(sampling_xy, "points",
                    atts = data.frame(dn = sampling_coord[, col_id]))
  terra::crs(p) <- terra::crs(r)

  # TODO (2026-09-04, existing limitation): vertex values pass through this
  # raster, so distinct subpixel positions in one cell may share a value.
  # Retained for this migration; replacing it needs a separate algorithm review.
  dummy_r <- terra::rasterize(p, r, field = "dn")

  mesh <- terra::delaunay(p)

  plane_coef <- lapply(1:length(mesh), function(i) {
    vertex_xy <- geom(mesh[i])[-1, c("x", "y")]
    z <- terra::extract(dummy_r, vertex_xy)
    pts <- cbind(vertex_xy, z = z[,])
    .plane_from_points(pts[1,], pts[2,], pts[3,])
  })
  plane_coef <- matrix(unlist(plane_coef), ncol = 3, byrow = TRUE)

  mesh <- cbind(mesh, plane_coef)
  names(mesh) <- c("alpha", "beta", "gamma")
  alpha <- rasterize(mesh, r, field = "alpha")
  beta <- rasterize(mesh, r, field = "beta")
  gamma <- rasterize(mesh, r, field = "gamma")

  x <- y <- r
  x[] <- xFromCell(r, 1:ncell(r))
  y[] <- yFromCell(r, 1:ncell(r))
  raster <- alpha* x + beta* y + gamma

  hull_xy <- sampling_xy[grDevices::chull(sampling_xy),]
  p <- terra::vect(hull_xy, "polygons", crs = crs(r))

  if (return == "raster") {
    return(raster)
  } else {
    return(list(raster = raster,
                mesh = mesh,
                convex_hull = p
                ))
  }
}






