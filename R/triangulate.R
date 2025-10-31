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
#' @note
#' The function internally uses [terra::delaunay()] to construct the
#' triangulated mesh.
#'
#' @inheritParams interpolate_planar
#'
#' @return Numeric [terra::SpatRaster-class] with one layer and the same geometry
#'   as `r`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' r <- caim$Blue
#'
#' bin <- binarize_by_region(r, ring_segmentation(z, 15), "thr_isodata") &
#'   select_sky_region(z, 0, 88)
#'
#' g <- sky_grid_segmentation(z, a, 10)
#' sky_points <- extract_sky_points(r, bin, g, dist_to_black = 3)
#' plot(bin)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#' sky_points <- extract_dn(r, sky_points)
#'
#' sky <- triangulate(sky_points, r, col_id = 3)
#' plot(sky)
#' }
triangulate <- function(sky_points, r, col_id = "dn") {
  if (!is.data.frame(sky_points)) {
    stop("`sky_points` must be a data frame.")
  }
  .assert_single_layer(r)
  handling_col_id <- c(
    tryCatch(.check_vector(col_id, "numeric", sign = "any"),
             error = function(e) FALSE),
    tryCatch(.check_vector(col_id, "character"),
             error = function(e) FALSE)
  )
  if (!any(handling_col_id)) {
    stop("`col_id` must be the name or position of the column in `sky_points` containing the values to triangulate.")
  }
  if (is.numeric(col_id)) col_id <- names(sky_points)[col_id]
  required_cols <- c("row", "col", col_id)

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

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  xy <- terra::xyFromCell(r, cells)

  p <- terra::vect(xy, "points", atts = data.frame(dn = sky_points[, col_id]))
  terra::crs(p) <- terra::crs(r)

  dummy_r <- terra::rasterize(p, r, field = "dn")

  mesh <- terra::delaunay(p)

  plane_coef <- lapply(1:length(mesh), function(i) {
    xy <- geom(mesh[i])[-1, c("x", "y")]
    z <- terra::extract(dummy_r, xy)
    pts <- cbind(xy, z = z[,])
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
  alpha* x + beta* y + gamma
}






