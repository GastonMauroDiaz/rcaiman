#' Fisheye to pano
#'
#' @param r
#' @param g
#' @param fun
 #' @noRd
setGeneric("fisheye_to_pano",
           function(r, z, a, fun = mean)
             standardGeneric("fisheye_to_pano"))

#' @rdname fisheye_to_pano
setMethod("fisheye_to_pano",
          signature(r = "RasterStackBrick"),
          function (r, z, a, fun)
          {
            layer_names <- names(r)
            pano <- Map(function(r) fisheye_to_pano(r, z, a, fun), as.list(r))
            pano <- stack(pano)
            names(pano) <- layer_names
            pano
          }
)

#' @rdname fisheye_to_pano
setMethod("fisheye_to_pano",
          signature(r = "RasterLayer"),
          function (r, z, a, fun)
          {
            g <- sky_grid_segmentation(z, a, 1)
            blue <- extract_feature(r, g, fun, return_raster = FALSE)
            xy <- .decode_label(as.numeric(names(blue)))
            r <- matrix(NA, ncol = max(xy$sector_ID), nrow = max(xy$rings_ID))
            r <- raster(r)
            extent(r) <- extent(0, ncol(r), 0, nrow(r))
            cells <- cellFromXY(r, as.matrix(xy) - 0.5)
            r[cells] <- blue
            flip(r, "y") #because nadir is 0 and for raster 0 is the bottom
          }
)
