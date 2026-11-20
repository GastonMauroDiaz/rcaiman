#' Extract a field-of-view contour
#'
#' @description
#' Extract the field-of-view contour from an internal geometry as an ordered
#' table of directions.
#'
#' @details
#' A pixel belongs to the border when both of its angles are present and at
#' least one of its four neighbors (above, below, left, or right) has a missing
#' angle or lies outside the image. The contour follows these border pixels
#' clockwise in the image, starting at the topmost, leftmost border pixel.
#' Ordering uses the image grid, not the numerical order of the angles.
#'
#' Each direction is taken directly from the center of a border pixel. The
#' contour therefore depends on the resolution of the input grid; it does not
#' describe the outer edges or the angular footprints of pixels. No directions
#' are interpolated, and no convex hull or simplification is applied.
#'
#' The first direction is repeated at the end to close the contour. Consecutive
#' occurrences of the same pixel are collapsed. A narrow connection or a
#' one-pixel-wide region can require revisiting a pixel, so the result is not
#' necessarily a simple polygon. An isolated valid pixel returns one direction.
#'
#' This function supports a single contour. Separate regions and holes produce
#' an error. Pixels touching only at a corner are treated as separate regions.
#' Masking the input changes the region whose contour is extracted.
#'
#' @param g numeric [terra::SpatRaster-class] with layers
#'   `off_axis` and `internal_azimuth`, in that order, in degrees, as returned
#'   by [set_internal_geometry()]. A missing angle excludes the pixel.
#'
#' @return data frame with columns `off_axis` and `internal_azimuth` in radians,
#'   ordered along the contour. If there are no valid directions, the table
#'   has zero rows.
#'
#' @seealso [set_internal_geometry()], [add_directions()],
#'   [intersect_fov()], [clip_fov()]
#' @export
#' @examples
#' geometry <- list(model = "radial_projection", parameters = equidistant_coef,
#'   center_of_vision = c(20, 15), normalization_radius = 18,
#'   normalization_angle = 90)
#' g <- set_internal_geometry(c(40, 30), geometry)
#' fov <- extract_fov(g)
#' head(fov)
extract_fov <- function(g) {
  angles <- .internal_coverage_angles(g, "g")
  cells <- .field_of_view_border(stats::complete.cases(angles),
                                 nrow(g), ncol(g))
  as.data.frame(angles[cells, , drop = FALSE] * pi / 180)
}

.field_of_view_border <- function(valid, height, width) {
  if (!any(valid)) return(integer())
  m <- matrix(valid, height, width, byrow = TRUE)
  padded <- matrix(FALSE, height + 2L, width + 2L)
  rows <- seq_len(height)
  cols <- seq_len(width)
  padded[rows + 1L, cols + 1L] <- m

  # Oriented grid edges keep the valid region on their right. Their owning
  # pixels supply the output directions; corner positions only order the walk.
  neighbors <- list(padded[rows, cols + 1L, drop = FALSE],
                    padded[rows + 1L, cols + 2L, drop = FALSE],
                    padded[rows + 2L, cols + 1L, drop = FALSE],
                    padded[rows + 1L, cols, drop = FALSE])
  edges <- lapply(0:3, function(direction) {
    rc <- which(m & !neighbors[[direction + 1L]], arr.ind = TRUE)
    x <- rc[, 2] - 1 + as.integer(direction %in% c(1L, 2L))
    y <- rc[, 1] - 1 + as.integer(direction %in% c(2L, 3L))
    from <- y * (width + 1) + x
    cbind(cell = (rc[, 1] - 1) * width + rc[, 2],
          direction = rep(direction, nrow(rc)), from = from,
          to = from + c(1, width + 1, -1, -width - 1)[direction + 1L])
  })
  edges <- do.call(rbind, edges)
  edges <- edges[order(edges[, "cell"], edges[, "direction"]), , drop = FALSE]
  key <- edges[, "from"] * 4 + edges[, "direction"]
  successor <- rep(NA_integer_, nrow(edges))
  # At diagonal contacts, turn right to stay on the same four-connected region.
  for (turn in c(1L, 0L, 3L)) {
    pending <- which(is.na(successor))
    successor[pending] <- match(edges[pending, "to"] * 4 +
      (edges[pending, "direction"] + turn) %% 4, key)
  }

  visited <- rep(FALSE, nrow(edges))
  path <- integer(nrow(edges))
  current <- 1L
  size <- 0L
  repeat {
    size <- size + 1L
    visited[current] <- TRUE
    path[size] <- edges[current, "cell"]
    current <- successor[current]
    if (current == 1L) break
  }
  if (!all(visited)) {
    stop("`g` contains more than one contour. ",
         "Supply a single region without holes.", call. = FALSE)
  }
  path <- path[seq_len(size)]
  path <- path[c(TRUE, diff(path) != 0)]
  if (length(path) > 1L && tail(path, 1L) != path[1L]) path <- c(path, path[1L])
  path
}
