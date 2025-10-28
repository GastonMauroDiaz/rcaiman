triangulate <- function(sky_points, r, col_id = "dn") {
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

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  xy <- terra::xyFromCell(r, cells)
  p <- terra::vect(xy, "points", atts = data.frame(dn = sky_points[, col_id]))
  terra::crs(p) <- terra::crs(r)

  dummy_r <- terra::rasterize(p, r, field = "dn")
  mesh <- terra::delaunay(p)

  l <- list()
  for (i in 1:length(mesh)) {
    xy <- geom(mesh[i])[-1, c("x", "y")]
    dn <- terra::extract(dummy_r, xy)
    rc  <- terra::rowColFromCell(dummy_r, cellFromXY(dummy_r, xy))
    dummy_sky_points <- data.frame(row = rc[,1], col = rc[,2], dn = dn[,])
    l[[i]] <- fit_trend_surface(dummy_sky_points, r, np = 1,
                                extrapolate = FALSE)$raster
  }

  stack <- rast(l)
  terra::cover(stack)
}






