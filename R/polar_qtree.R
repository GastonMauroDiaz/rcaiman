#' Do quad-tree segmentation in the polar space
#'
#' The quad-tree segmentation algorithm is a top-down process that makes
#' recursive divisions in four equal parts until a condition is satisfied and
#' stops locally. The usual implementation of the quad-tree algorithm is
#' based on the raster structure and this is why the result are squares of
#' different sizes. This method implements the quad-tree segmentation in a polar
#' space, so the segments are shaped like windshields, though some of them will
#' look elongated in height. The pattern is two opposite and converging straight
#' sides and two opposite and parallel curvy sides.
#'
#' The algorithm splits segments of 30 degrees resolution into four sub-segments
#' and calculates the standard deviation of the pixels from `r` delimited
#' by each of those segments. The splitting process stops locally if the sum of
#' the standard deviation of the sub-segments minus the standard deviation of
#' the parent segment (named *delta*) is less or equal than the
#' `scale_parameter`. If `r` has more than one layer, *delta* is
#' calculated separately and *delta* mean is used to evaluate the stopping
#' condition.
#'
#' @param r [SpatRaster-class].
#' @inheritParams sky_grid_segmentation
#' @param scale_parameter Numeric vector of length one. Quad-tree is a top-down
#'   method. This parameter controls the stopping condition. Therefore, it
#'   allows controlling the size of the resulting segments. Ultimately, segments
#'   sizes will depend on both this parameter and the heterogeneity of `r`.
#'
#' @return A single layer image of the class [SpatRaster-class] with
#'   integer values.
#'
#' @export polar_qtree
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize_minmax()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' seg <- polar_qtree(caim, z, a)
#' plot(seg)
#' plot(extract_feature(caim$Blue, seg))
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
  ges <- Map(sky_grid_segmentation, z, a, angle.wds)

  .calc_delta_single_layer <- function(r) {
    if (any(is.na(r)[] %>% as.logical())) {
      .sd <- function(x) {
        x <- sd(x, na.rm = TRUE)
        if (is.na(x)) x <- 0
        x
      }
      sd_now <- Map(function(g) extract_feature(r, g, .sd,
                                                return_raster = FALSE),
                    ges)
    } else {
      sd_now <- Map(function(g) extract_feature(r, g, sd,
                                                return_raster = FALSE),
                    ges)
    }
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
  seg <- max(seg)
  names(seg) <- "Polar quad-tree"
  seg
}
