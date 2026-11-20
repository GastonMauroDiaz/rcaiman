#' Create a convex hull mask from image coordinates
#'
#' @description
#' Delimit the convex hull of sample positions on an image grid.
#'
#' @details
#' The hull is computed from `x` and `y_down`, independently of sample values
#' or a fitted model. Duplicate positions are allowed; at least three distinct,
#' non-collinear positions are needed to define a polygon.
#'
#' Rasterization uses [terra::rasterize()] with `touches = FALSE`, selecting
#' cells by their centers rather than by any overlap with the hull. The grid
#' clips the polygon; samples outside the grid still contribute to its shape.
#' Template values, including `NA`, are ignored, so the mask can be combined
#' explicitly with any other selection.
#'
#' @param sampling_coord `data.frame` with finite numeric columns `x` and
#'   `y_down`, in the coordinate frame of `r`. Additional columns are ignored.
#' @param r [terra::SpatRaster-class] with one layer and the unit-resolution
#'   grid established by [set_rcaiman_geometry()]. Only its grid is used.
#'
#' @return logical [terra::SpatRaster-class] with one layer named
#'   `convex_hull`, `TRUE` for cells inside the rasterized hull and `FALSE`
#'   elsewhere. The mask contains no `NA` values.
#'
#' @seealso [fit_trend_surface()], [trend_surface_image()]
#' @export
#' @examples
#' r <- set_rcaiman_geometry(terra::rast(nrows = 8, ncols = 10))
#' samples <- data.frame(x = c(1, 9, 2), y_down = c(1, 2, 7))
#' m <- convex_hull_mask(samples, r)
#' plot(m)
#' points(samples$x, y_from_y_down(samples$y_down, nrow(r)), col = 2)
convex_hull_mask <- function(sampling_coord, r) {
  .check_columns(sampling_coord, required_columns = c("x", "y_down"))
  .check_vector(sampling_coord$x, "numeric")
  .check_vector(sampling_coord$y_down, "numeric")
  .check_image_template(r)
  positions <- unique(sampling_coord[c("x", "y_down")])
  hull <- if (nrow(positions) >= 3) grDevices::chull(positions) else integer()
  if (length(hull) < 3) {
    stop("A convex hull requires at least three non-collinear positions.",
         call. = FALSE)
  }
  xy <- cbind(positions$x[hull], terra::nrow(r) - positions$y_down[hull])
  polygon <- terra::vect(xy, type = "polygons", crs = terra::crs(r))
  m <- terra::rasterize(polygon, r, field = 1, background = 0,
                         touches = FALSE) == 1
  names(m) <- "convex_hull"
  m
}
