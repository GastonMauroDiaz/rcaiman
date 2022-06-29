#' Quad-tree segmentation in a polar space.
#'
#' The quad-tree segmentation algorithm is a top-down process that makes
#' recursive divisions in four equal parts until a condition is satisfied and
#' then stops locally. The usual implementation of the quad-tree algorithm is
#' based on the raster structure and this is why the result are squares of
#' different sizes. This method implements the quad-tree segmentation in a polar
#' space, so the segments are shaped like windshields (please see
#' \code{\link{sky_grid_segmentation}}).
#'
#' The algorithm splits segments of 30 degrees resolution into four sub-segments
#' and calculates the standard deviation of the pixels from \code{r} delimited
#' by each of those segments. The splitting process stops locally if the sum of
#' the standard deviation of the sub-segments minus the standard deviation of
#' the parent segment (named \emph{delta}) is less or equal than the
#' \code{scale_parameter}. If \code{r} has more than one layer, \emph{delta} is
#' calculated separately and \emph{delta} mean is used to evaluate the stopping
#' condition.
#'
#' @param r \linkS4class{SpatRaster}.
#' @inheritParams ootb_mblt
#' @param scale_parameter Numeric vector of length one. Quad-tree is a top-down
#'   method. This parameter controls the stopping condition. Therefore, it
#'   allows the user to controls the size of the resulting segments. See Details
#'   for more information.
#'
#' @export polar_qtree
#'
#' @family Segmentation functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' caim <- normalize(caim, 0, 255)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' seg <- polar_qtree(caim, z, a)
#' plot(seg)
#' tmp <- extract_feature(caim$Blue, seg)
#' plot(tmp)
#' }
polar_qtree <- function(r, z, a,
                        scale_parameter = 0.2) {
  stopifnot(class(r) == "SpatRaster")
  .is_single_layer_raster(z)
  .is_single_layer_raster(a)
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
  terra::compareGeom(z, a)

  # mnSize <- 1000
  angle.wds <- c(15, 7.5, 3.75, 1.875)

  zes <- Map(function(i) z, seq_along(angle.wds))
  aes <- Map(function(i) a, seq_along(angle.wds))
  ges <- Map(sky_grid_segmentation, zes, aes, angle.wds)

  .calc_delta_single_layer <- function(r) {
    sd_now <- Map(function(g) extract_feature(r, g, sd,
                                              return_raster = FALSE),
                  ges)
    indices_if_split <- Map(function(i) extract_feature(ges[[i-1]], ges[[i]],
                                                        max,
                                                        return_raster = FALSE),
                            2:length(ges))
    sd_if_split <- Map(function(i) tapply(sd_now[[i+1]],
                                          indices_if_split[[i]], sum),
                       1:(length(ges)-1))
    delta <- Map(function(i) sd_if_split[[i]] - sd_now[[i]],
                 seq_along(sd_if_split))
    delta
  }

  if (terra::nlyr(r) > 1) {
    delta <- Map(.calc_delta_single_layer, as.list(r))
    delta <- Map(function(i) {
                  x <- Map(function(j) delta[[j]][[i]], seq_along(delta))
                  apply(as.data.frame(x), 1, mean)
               }, seq_along(delta[[1]]))
  } else {
    delta <- .calc_delta_single_layer(r)
  }
  it_should_be_splited <- Map(function(x) x > scale_parameter, delta)
  it_should_be_splited <- Map(function(i) {
    terra::subst(ges[[i]],
                 names(delta[[i]]) %>%
                   as.numeric(),
                 it_should_be_splited[[i]])
  }, seq_along(it_should_be_splited))
  ges <- ges[-1]
  ges <- Map(function(i) {
    r <- ges[[i]]
    r[r == 0] <- NA
    r + i*10000000
  }, seq_along(ges))
  ges <- terra::rast(ges)
  it_should_be_splited <- terra::rast(it_should_be_splited)
  seg <- ges * it_should_be_splited
  g <- sky_grid_segmentation(z, a, angle.wds[1])
  seg <- c(g, seg)
  max(seg)
}
