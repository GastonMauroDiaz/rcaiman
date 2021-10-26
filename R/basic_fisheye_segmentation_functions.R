
#' Rings segmentation
#'
#' Segmenting an hemispherical view by slicing the zenith angle from \code{0}
#' to \code{90º} in equals intervals.
#'
#' @inheritParams azimuth_image
#' @param angle_width Length-one numeric vector. Angle in degrees able to
#'   divide the angle range into a whole number of segments.
#' @param return_angle Length-one logical vector. If it is \code{FALSE}, all the
#'   pixels that belong to a segment are labeled with an ID number.
#'   Otherwise, the angle mean of the segment is assigned to the pixels.
#'
#' @return \linkS4class{RasterLayer} with segments shaped like concentric rings.
#' @export
#'
#' @family basic fisheye segmentation functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' rings <- rings_segmentation(z, 15)
#' plot(rings == 1)
rings_segmentation <- function(z, angle_width, return_angle = FALSE) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(return_angle) == "logical")
  stopifnot(length(angle_width) == 1)

  if (!.is_whole(90 / angle_width)) {
    stop(
      paste("angle_width should divide the,",
            "0 to 90 range into a whole number of segments.")
    )
  }

  intervals <- seq(0, 90, angle_width)
  c1 <- intervals[1:(length(intervals) - 1)]
  c2 <- intervals[2:length(intervals)]

  if (return_angle) {
    c3 <- (c1 + c2) / 2
  } else {
    c3 <- 1:(length(intervals) - 1)
  }
  rcl <- matrix(c(c1, c2, c3), ncol = 3)
  reclassify(z, rcl)
}



#' Sectors segmentation
#'
#' Segmenting an hemispherical view by slicing the azimuth angle from \code{0}
#' to \code{360º} in equals intervals.
#'
#' @param a \code{\linkS4class{RasterLayer}} built with
#'   \code{\link{azimuth_image}}.
#' @inheritParams rings_segmentation
#'
#' @return \linkS4class{RasterLayer} with segments shaped like pizza slices.
#' @export
#'
#' @family basic fisheye segmentation functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' a <- azimuth_image(z)
#' sectors <- sectors_segmentation(a, 15)
#' plot(sectors == 1)
sectors_segmentation <- function(a, angle_width, return_angle = FALSE) {
  stopifnot(class(a) == "RasterLayer")
  stopifnot(class(return_angle) == "logical")
  stopifnot(length(angle_width) == 1)

  if (!.is_whole(360 / angle_width)) {
    stop(
      paste("angle_width should divide the 0 to 360",
            "range into a whole number of segments.")
    )
  }

  intervals <- seq(0, 360, angle_width)
  c1 <- intervals[1:(length(intervals) - 1)]
  c2 <- intervals[2:length(intervals)]
  if (return_angle) {
    c3 <- (c1 + c2) / 2
  } else {
    c3 <- 1:(length(intervals) - 1)
  }

  rcl <- matrix(c(c1, c2, c3), ncol = 3)
  reclassify(a, rcl)
}



#' Sky grid segmentation
#'
#' Segmenting the hemisphere view into segments of equal angular resolution for
#' both zenith and azimuth angles.
#'
#' Intersecting rings with sectors makes a grid in which each segment is a
#' portion of the hemisphere. Each pixel of the grid is labeled with an ID that
#' codify both ring and sector ID. For example, a grid with a regular interval
#' of 1 degree has segment from \code{1001} to \code{360090}. This numbers are
#' calculated with: \code{sectorID * 1000 + ringsID}, where \code{sectorID} is
#' the ID number of the sector and \code{ringsID} is the ID number of the ring.
#'
#'
#' @inheritParams azimuth_image
#' @inheritParams sectors_segmentation
#' @param angle_width One-length numeric vector. It should be \code{30, 15, 10,
#'   7.5, 6, 5, 3.75, 3, 2.5, 1.875, 1} or \code{0.5} degrees. This constrains
#'   is rooted in the requirement of a value able to divide both the \code{0} to
#'   \code{360} and \code{0} to \code{90} ranges into a whole number of
#'   segments.
#' @param sequential One-length logical vector. If it is \code{TRUE} the
#'   segments are labeled with sequential numbers. By default (\code{FALSE}),
#'   labeling numbers are not sequential (see Details).
#'
#' @return \linkS4class{RasterLayer} with segments shaped like windshields,
#'   although some of them will look elongated in height. The pattern is two
#'   opposite and converging straight sides and two opposite and parallel
#'   curvy sides.
#' @export
#'
#' @family basic fisheye segmentation functions
#'
#' @examples
#' z <- zenith_image(1490, lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 15)
#' plot(g == 36009)
#' g <- sky_grid_segmentation(z, a, 15, sequential = TRUE)
#' plot(g, col = sample(rainbow(length(unique(g)))))
sky_grid_segmentation <- function(z, a, angle_width, sequential = FALSE) {
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(.get_max(z) < 90)

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

  g
}
