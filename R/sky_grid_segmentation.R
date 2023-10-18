#' Sky grid segmentation
#'
#' Segmenting the hemisphere view into segments of equal angular resolution for
#' both zenith and azimuth angles.
#'
#' Intersecting rings with sectors makes a grid in which each cell is a
#' portion of the hemisphere. Each pixel of the grid is labeled with an ID that
#' codify both ring and sector IDs. For example, a grid with a regular interval
#' of one degree has segment from `1001` to `360090`. This numbers are
#' calculated with: \eqn{sectorID \times 1000 + ringID}, where \eqn{sectorID} is
#' the ID number of the sector and \eqn{ringID} is the ID number of the ring.
#'
#'
#' @inheritParams ootb_mblt
#' @param angle_width Numeric vector of length one. It should be `30, 15,
#'   10, 7.5, 6, 5, 3.75, 3, 2.5, 1.875, 1` or `0.5` degrees. This
#'   constrain is rooted in the requirement of a value able to divide both the
#'   `0` to `360` and `0` to `90` ranges into a whole number
#'   of segments.
#' @param sequential Logical vector of length one. If it is `TRUE`, the
#'   segments are labeled with sequential numbers. By default (`FALSE`),
#'   labeling numbers are not sequential (see Details).
#'
#' @return An object from the class [SpatRaster-class] with segments
#'   shaped like windshields, though some of them will look elongated in
#'   height. The pattern is two opposite and converging straight sides and two
#'   opposite and parallel curvy sides.
#' @export
#'
#' @family Segmentation Functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 15)
#' plot(g == 24005)
#' \dontrun{
#' g <- sky_grid_segmentation(z, a, 15, sequential = TRUE)
#' col <- terra::unique(g) %>% nrow() %>% rainbow() %>% sample()
#' plot(g, col = col)
#' }
sky_grid_segmentation <- function(z, a, angle_width, sequential = FALSE) {
  stopifnot(length(sequential) == 1)
  stopifnot(class(sequential) == "logical")
  stopifnot(length(angle_width) == 1)

  if (!max(angle_width == c(30,  15,    10,   7.5,
                            6,   5,     3.75, 3,
                            2.5, 1.875, 1,    0.5))) {
    stop(
      paste0("angle_width should be 30, 15, 10, 7.5,",
             "6, 5, 3.75, 3, 2.5, 1.875, 1 or 0.5 degrees.")
    )
  }

  fun <- function(s, r) s * 1000 + r
  g <- fun(sectors_segmentation(a, angle_width),
           rings_segmentation(z, angle_width))

  if (sequential) {
    from <- unique(terra::values(g))
    to <- 1:length(from)
    g <- terra::subst(g, from, to)
  }
  names(g) <- "Sky grid"
  g
}
