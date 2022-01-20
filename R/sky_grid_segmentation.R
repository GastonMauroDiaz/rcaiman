#' Sky grid segmentation
#'
#' Segmenting the hemisphere view into segments of equal angular resolution for
#' both zenith and azimuth angles.
#'
#' Intersecting rings with sectors makes a grid in which each segment is a
#' portion of the hemisphere. Each pixel of the grid is labeled with an ID that
#' codify both ring and sector IDs. For example, a grid with a regular interval
#' of one degree has segment from \code{1001} to \code{360090}. This numbers are
#' calculated with: \code{sectorID * 1000 + ringsID}, where \code{sectorID} is
#' the ID number of the sector and \code{ringsID} is the ID number of the ring.
#'
#'
#' @inheritParams azimuth_image
#' @inheritParams sectors_segmentation
#' @param angle_width Numeric vector of length one. It should be \code{30, 15,
#'   10, 7.5, 6, 5, 3.75, 3, 2.5, 1.875, 1} or \code{0.5} degrees. This
#'   constrain is rooted in the requirement of a value able to divide both the
#'   \code{0} to \code{360} and \code{0} to \code{90} ranges into a whole number
#'   of segments.
#' @param sequential Logical vector of length one. If it is \code{TRUE}, the
#'   segments are labeled with sequential numbers. By default (\code{FALSE}),
#'   labeling numbers are not sequential (see Details).
#'
#' @return An object from the class \linkS4class{RasterLayer} with segments
#'   shaped like windshields, although some of them will look elongated in
#'   height. The pattern is two opposite and converging straight sides and two
#'   opposite and parallel curvy sides.
#' @export
#'
#' @family Segmentation functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 15)
#' plot(g == 24005)
#' \dontrun{
#' g <- sky_grid_segmentation(z, a, 15, sequential = TRUE)
#' plot(g, col = sample(rainbow(length(raster::unique(g)))))
#' }
sky_grid_segmentation <- function(z, a, angle_width, sequential = FALSE) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(.get_max(z) <= 90)

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

  g <- overlay(
    sectors_segmentation(a, angle_width),
    rings_segmentation(z, angle_width),
    fun = fun
  )

  if (sequential) {
    df <- levels(as.factor(g))[[1]]
    df <- cbind(df, 1:nrow(df))
    g <- raster::subs(g, df)
  }

  as.factor(g)
}
