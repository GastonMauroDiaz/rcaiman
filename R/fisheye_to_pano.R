#' Fisheye to pano
#'
#' @param r
#' @param g
#' @param m
#' @param fill_na
#' @param fun
 #' @noRd
setGeneric("fisheye_to_pano",
           function(r, g, m = NULL, fill_na = TRUE, fun = mean)
             standardGeneric("fisheye_to_pano"))

#' @rdname fisheye_to_pano
setMethod("fisheye_to_pano",
          signature(r = "RasterStackBrick"),
          function (r, g, m, fill_na, fun)
          {
            layer_names <- names(r)
            pano <- Map(function(r) {
                                      fisheye_to_pano(r, g, m, fill_na, fun)
                                    }, as.list(r))
            pano <- stack(pano)
            names(pano) <- layer_names
            pano
          }
)

#' @rdname fisheye_to_pano
setMethod("fisheye_to_pano",
          signature(r = "RasterLayer"),
          function (r, g, m, fill_na, fun)
          {
            if(!is.null(m)) r[!m] <- NA
            blue <- extractFeatures(r, g, fun)
            xy <- .decode_label(as.numeric(names(blue)))
            r <- matrix(NA, ncol = max(xy$sector_ID), nrow = max(xy$rings_ID))
            r <- raster(r)
            extent(r) <- extent(0, ncol(r), 0, nrow(r))

            cells <- cellFromXY(r, as.matrix(xy)-0.5)
            r[cells] <- blue
            r <- flip(r, "y") #because nadir is 0 and for raster 0 is the bottom

            if(!is.null(m)) r <- trim(r)

            if (fill_na) {
              fun <- function(x) {
                x <- focal(x, matrix(1, 5, 5), mean, NAonly = TRUE, na.rm = TRUE)
                no_of_loops <- 0
                unlock <- TRUE
                while (unlock) {
                  no_of_loops <- no_of_loops + 1
                  x <- focal(x, matrix(1, 5, 5), mean, NAonly = TRUE, na.rm = TRUE)
                  unlock <- any(is.na(x)[])
                  if (unlock) unlock <- no_of_loops < 100
                }
                x
              }

              r <- fun(r)
            }
            r
          }
)
