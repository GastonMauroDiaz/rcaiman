#' Quad-tree segmentation in a polar space.
#'
#' The quad-tree segmentation algorithm is a top-down process that makes
#' recursive divisions in four equal parts until a condition is satisfied and
#' then stops locally. The usual implementation of the quad-tree algorithm is
#' based on the raster structure and this is why the result are squares of
#' different sizes. This method implements the quad-tree segmentation in a polar
#' space.
#'
#' Argument division could be 1, 2, 3 or 4 and it controls segment resolution.
#'
#' value / resolution (degrees)
#'
#' 1 / 30, 15, 7.5, 3.75 or 1.875.
#'
#' 2 / 15, 7.5, 3.75 or 1.875.
#'
#' 3 / 15, 7.5 or 3.75.
#'
#' 4 / 7.5 or 3.75.
#'
#' The algorithm starts with segments of a given resolution depending on
#' \code{divisions}. Next, it selects one segment and splits it into four
#' segments of equal angular resolution. Then, it uses the standard deviation of
#' \code{x} as homogeneous criterion. To that end, it calculates the standard
#' deviation for the entire segment and for each four segments. To stop the
#' process locally, the algorithm evaluates if the sum of the standard deviation
#' of the subsegments minus the standard deviation of the segment (delta) is
#' less or equal than the \code{scaleParameter}. If x is multilayer delta is
#' calculated separately and delta mean is used to evaluate the stopping
#' condition.
#'
#' @param x \code{\linkS4class{Raster}}. The raster to be processed. Single or
#'   multi-layer.
#' @param z Should match the lens geometry of
#'   the picture linked with \code{x}.
#' @param a raster
#' @param scaleParameter one-length numeric. It is part of the stopping
#'   condition (see Detalis).
#' @param divisions numeric. See details
#' @param mnSize numeric.  This is an additional stopping criterion to avoid
#'   errors in segments closed to zenith.
#' @param parallel logical. Go parallel.
#' @param freeThreads numeric. The number of threads that remain free for other
#'   tasks. Defaul \code{2} avoid the computer comes to a standstill until the R
#'   task has finished. Using \code{1} could be a good idea. Use \code{0} at
#'   your own risk. If you run out of memory, increase \code{freeThreads}.
#' @param memoryUseFactor numeric. Increase to use less memory when
#'   \code{parallel} is set to TRUE.
#'
#' @export polar_qtree
#' @noRd
polar_qtree <- function(x, z, a,
                        scaleParameter = 2,
                        divisions = 1,
                        mnSize = 1000,
                        parallel = TRUE,
                        freeThreads = 0) {
  print("please wait... performing quad-tree segmentation.")

  stopifnot(length(divisions) == 1)
  stopifnot(any(
    divisions == 1, divisions == 2,
    divisions == 3, divisions == 4
  ))

  raster::compareRaster(x, z)
  raster::compareRaster(z, a)

  cropFun <- function(x, cropped) {
    x <- raster::crop(x, cropped) * cropped
    x[cropped == 0] <- NA
    x
  }

  angle.wds <- list(
    c(15, 7.5, 3.75, 1.875),
    c(7.5, 3.75, 1.875),
    c(7.5, 3.75),
    c(3.75)
  )
  angle.wds <- angle.wds[[divisions]]

  g <- sky_grid_segmentation(z, a, max(angle.wds) * 2)

  if (parallel) {
    # Calculate the number of threads
    no_threads <- parallel::detectCores() - freeThreads

    rings <- rings_segmentation(z, 30)
    sectors <- sectors_segmentation(a, 60)
    polarGrid <- as.factor(sectors * 1000 + rings)
    segments <- levels(polarGrid)[[1]][, 1]
    ges <- Map(function(x) as.factor(g[polarGrid == x, drop = FALSE]), segments)

    # go parallel
    ## Initiate cluster
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl, c("x", "z", "a", "scaleParameter", "cropFun", "angle.wds"), environment())
    ges <- parallel::parLapply(cl, ges, function(g) {
      segments <- levels(g)[[1]][, 1]
      segments <- Map(function(x) x, segments)
      goRecursively <- function(x, z, a, segment, scaleParameter) {
        cropped <- g[g == segment, drop = FALSE]
        cropped <- !is.na(cropped)

        x <- cropFun(x, cropped)
        z <- cropFun(z, cropped)
        a <- cropFun(a, cropped)

        angleResolution <- (.get_max(z) - .get_min(z) + .get_max(a) - .get_min(a)) / 2

        if (angleResolution > min(angle.wds)) {
          angle.wd <- angle.wds[which.min(abs(angle.wds - rep(
            angleResolution / 2,
            length(angle.wds)
          )))]

          g2 <- sky_grid_segmentation(z, a, angle.wd)

          segments2 <- levels(g2)[[1]][, 1]

          if (length(segments2) > 4) stop("Please, report error code qt01")
          if (length(segments2) < 4) {
            delta <- scaleParameter - 1
          } else {
            hfun <- function(x) stats::sd(x, na.rm = TRUE)

            if (any(class(x)[1] == "RasterStack", class(x)[1] == "RasterBrick")) {
              delta <- c()
              for (l in 1:nlayers(x)) {
                sdIfSplit <- sum(extract_feature(raster::subset(x, l), g2, hfun, return_raster = FALSE))
                sdNow <- hfun(values(raster::subset(x, l)))
                delta[l] <- sdIfSplit - sdNow
              }

              delta <- sum(delta) / nlayers(x)
            } else {
              sdIfSplit <- sum(extract_feature(x, g2, hfun, return_raster = FALSE))
              sdNow <- hfun(values(x))
              delta <- sdIfSplit - sdNow
            }
          }

          if (delta > scaleParameter) {
            if (all(nrow(levels(g2)[[1]]) == 4, raster::ncell(g2) > mnSize)) {
              segments <- 1:4 + .get_max(g)
              df <- data.frame(levels(g2), 1:4)
              g2 <- raster::subs(g2, df)
              index <- !is.na(values(g2))
              g[g == segment] <<- values(g2)[index] + .get_max(g)

              segments <- Map(function(x) x, segments)
              baz <- function(segment) goRecursively(x, z, a, segment, scaleParameter)

              lapply(segments, baz)
            }
          }
        }
      }
      lapply(segments, function(segment) goRecursively(x, z, a, segment, scaleParameter))
      g
    })

    ## finish
    parallel::stopCluster(cl)

    # to avoid segments with same id
    maxOfGes <- unlist(Map(.get_max, ges))

    for (i in 2:length(maxOfGes)) {
      maxOfGes[i] <- maxOfGes[i] + maxOfGes[i - 1]
    }
    maxOfGes[1] <- 0
    ges <- Map(function(x, y) x + y, ges, maxOfGes)

    # merge
    g[] <- NA

    for (i in 1:(length(ges))) {
      g <- cover(g, extend(ges[[i]], g))
    }
  } else {
    # I repeat myself because environment
    goRecursively <- function(x, z, a, segment, scaleParameter) {
      cropped <- g[g == segment, drop = FALSE]
      cropped <- !is.na(cropped)

      x <- cropFun(x, cropped)
      z <- cropFun(z, cropped)
      a <- cropFun(a, cropped)

      angleResolution <- (.get_max(z) - .get_min(z) + .get_max(a) - .get_min(a)) / 2

      if (angleResolution > min(angle.wds)) {
        angle.wd <- angle.wds[which.min(abs(angle.wds - rep(
          angleResolution / 2,
          length(angle.wds)
        )))]

        g2 <- sky_grid_segmentation(z, a, angle.wd)

        segments2 <- levels(g2)[[1]][, 1]

        if (length(segments2) > 4) stop("Please, report error code qt01")
        if (length(segments2) < 4) {
          delta <- scaleParameter - 1
        } else {
          hfun <- function(x) stats::sd(x, na.rm = TRUE)

          if (any(class(x)[1] == "RasterStack", class(x)[1] == "RasterBrick")) {
            delta <- c()
            for (l in 1:nlayers(x)) {
              sdIfSplit <- sum(extract_feature(raster::subset(x, l), g2, hfun, return_raster = FALSE))
              sdNow <- hfun(values(raster::subset(x, l)))
              delta[l] <- sdIfSplit - sdNow
            }

            delta <- sum(delta) / nlayers(x)
          } else {
            sdIfSplit <- sum(extract_feature(x, g2, hfun, return_raster = FALSE))
            sdNow <- hfun(values(x))
            delta <- sdIfSplit - sdNow
          }
        }

        if (delta > scaleParameter) {
          if (all(nrow(levels(g2)[[1]]) == 4, raster::ncell(g2) > mnSize)) {
            segments <- 1:4 + .get_max(g)
            df <- data.frame(levels(g2), 1:4)
            g2 <- raster::subs(g2, df)
            index <- !is.na(values(g2))
            g[g == segment] <<- values(g2)[index] + .get_max(g)

            segments <- Map(function(x) x, segments)
            baz <- function(segment) goRecursively(x, z, a, segment, scaleParameter)

            lapply(segments, baz)
          }
        }
      }
    }
    segments <- levels(g)[[1]][, 1]
    pb <- raster::pbCreate(length(segments), "text")
    for (i in 1:length(segments)) {
      raster::pbStep(pb, i)
      goRecursively(x, z, a, segments[i], scaleParameter)
    }
    raster::pbClose(pb)
  }
  g
}
