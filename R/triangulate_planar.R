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
#' the interpolated surface in a meaninful way.
#'
#' No consistency checks are performed to ensure that `sampling_points` and `r`
#' are geometrically compatible. Incorrect combinations may lead to invalid
#' outputs.
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
#'
#' sky <- triangulate_planar(sky_points, "Cyan", r)
#' plot(sky)
#' }
triangulate_planar <- function(sky_points, col_id = "dn", r, return = "raster") {
  if (!is.data.frame(sky_points)) {
    stop("`sky_points` must be a data frame.")
  }
  .assert_single_layer(r, allow_na = FALSE)
  required_cols <- .resolve_required_column(col_id, sky_points)
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
  raster <- alpha* x + beta* y + gamma

  xy <- xy[grDevices::chull(xy),]
  p <- terra::vect(xy, "polygons", crs = crs(r))

  if (return == "raster") {
    return(raster)
  } else {
    return(list(raster = raster,
                mesh = mesh,
                convex_hull = p
                ))
  }
}






