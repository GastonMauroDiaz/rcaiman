#' Do sky grid segmentation
#'
#' Segment the hemisphere view into segments of equal angular resolution for
#' both zenith and azimuth angles.
#'
#' Intersecting rings with sectors makes a grid in which each cell is a portion
#' of the hemisphere. Each pixel of the grid is labeled with an ID that codify
#' both ring and sector IDs. For example, a grid with a regular interval of one
#' degree has segment from `1001` to `360090`. This numbers are calculated with:
#' \eqn{sectorID \times 1000 + ringID}, where \eqn{sectorID} is the ID number of
#' the sector and \eqn{ringID} is the ID number of the ring.
#'
#'
#' @param z [SpatRaster-class] built with [zenith_image()].
#' @param a [SpatRaster-class] built with [azimuth_image()].
#' @param angle_width Numeric vector of length one. It should be `30, 15, 10,
#'   7.5, 6, 5, 3.75, 3, 2.5, 1.875, 1` or `0.5` degrees. This constrain is
#'   rooted in the requirement of a value able to divide both the `0` to `360`
#'   and `0` to `90` ranges into a whole number of segments.
#' @param first_ring_different Logical vector of length one. If it is `TRUE`,
#'   the first ring (the one containing the zenith) is not divided into cells.
#' @param sequential Logical vector of length one. If it is `TRUE`, the segments
#'   are labeled with sequential numbers. By default (`FALSE`), labeling numbers
#'   are not sequential (see Details).
#'
#' @return An object of the class [SpatRaster-class] with segments shaped like
#'   windshields, though some of them will look elongated in height. The pattern
#'   is two opposite and converging straight sides and two opposite and parallel
#'   curvy sides.
#' @export
#'
#' @examples
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' g <- sky_grid_segmentation(z, a, 15)
#' plot(g == 24005)
#' \dontrun{
#' # This code output the angles that can be used. For convenience, the
#' # column radians_denom can be used to provide `angle_width` as:
#' # 180/[value from radians_denom]
#' df <- data.frame(degrees = 90 / 1:180)
#'
#' deg_to_pi_expr <- function(deg) {
#'   frac <- MASS::fractions(deg / 180)
#'   strsplit(as.character(frac), "/")[[1]][2] %>% as.numeric()
#' }
#'
#' df$radians_denom <- sapply(df$degrees, function(deg) deg_to_pi_expr(deg))
#'
#' z <- zenith_image(10, lens())
#' a <- azimuth_image(z)
#' u <- c()
#' for (i in 1:nrow(df)) {
#'   u <- c(u, tryCatch(is((sky_grid_segmentation(z, a,
#'                             180/df$radians_denom[i])), "SpatRaster"),
#'                      error = function(e) FALSE))
#' }
#' df <- df[u,]
#' df
#'
#' # here an example
#' g <- sky_grid_segmentation(z, a, 180/14, sequential = TRUE)
#' col <- terra::unique(g) %>% nrow() %>% rainbow() %>% sample()
#' plot(g, col = col)
#'
#' g <- terra::focal(g, 3, sd)
#' g <- g != 0
#' plot(g)
#' plot(!g * caim$Blue)
#'
#' }
sky_grid_segmentation <- function(z, a, angle_width,
                                  first_ring_different = FALSE,
                                  sequential = FALSE) {
  stopifnot(length(sequential) == 1)
  stopifnot(class(sequential) == "logical")
  stopifnot(length(angle_width) == 1)

  # if (!max(angle_width == c(30,  15,    10,   7.5,
  #                           6,   5,     3.75, 3,
  #                           2.5, 1.875, 1,    0.5))) {
  #   stop(
  #     paste0("angle_width should be 30, 15, 10, 7.5, ",
  #            "6, 5, 3.75, 3, 2.5, 1.875, 1 or 0.5 degrees.")
  #   )
  # }

  fun <- function(s, r) s * 1000 + r
  g <- fun(sectors_segmentation(a, angle_width),
           rings_segmentation(z, angle_width))

  if (first_ring_different) {
    for (i in seq(1, 360/angle_width)*1000 + 1) {
      g[g == i] <- 1000
    }
  }

  if (sequential) {
    from <- unique(terra::values(g))
    to <- 1:length(from)
    g <- terra::subst(g, from, to)
  }
  names(g) <- paste0("Sky grid, ", angle_width, " degrees")
  g
}
