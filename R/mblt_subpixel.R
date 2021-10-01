.interpolate_akima <- function(x, linear = TRUE) {

  xy <- coordinates(x)
  z <- values(x)

  xy <- xy[!is.na(z),]
  z <- z[!is.na(z)]

  xo <- seq(xmin(x), xmax(x), length=ncol(x))
  yo <- seq(ymin(x), ymax(x), length=nrow(x))
  r <- akima::interp(xy[,1], xy[,2], z, xo=xo, yo=yo, linear)

  r <- raster(r)
  extent(r) <- extent(x)
  r <- extend(r, x)
  projection(r) <- projection(x)
  r
}

#' MBLT subpixel
#'
#'
#' @inheritParams fit_trens_surface_to_sky_dn
#' @param bin_plant \linkS4class(RasterLayer). A binarization of r that works
#'   as a mask for pure plant pixels.
#' @param bin_sky \linkS4class(RasterLayer). A binarization of r that works
#'   as a mask for pure sky pixels.
#' @inheritParams model_sky_dn
#'
#' @return \linkS4class(RasterLayer)
#'
mblt_subpixel <- function (r, m, bin_plant, bin_sky,
                           parallel = TRUE,
                           free_threads = 0) {

  if (!requireNamespace("akima", quietly = TRUE)) {
    stop(paste("Package \"akima\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  sample_size <-  30

  compareRaster(r, m)
  compareRaster(r, bin_plant)
  compareRaster(r, bin_sky)
  if (!is.null(fill)) compareRaster(r, fill)

  blue <- r * 255
  rm(r)

  size_range <- ncol(blue) * 0.02
  size_range <- round(size_range)
  if(!.is_even(size_range)) size_range <- size_range + 1

  size_step <- round(size_range * 0.2)
  if (!.is_even(size_step)) size_step <- size_step - 1

  bin_plant[!m] <- 0
  bin_sky[!m] <- 0

  sizes <- round(sqrt(sample_size*4))
  if (is.even(sizes)) sizes <- sizes + 1
  sizes <- seq(sizes, sizes + size_range, size_step)


  fun <- function(x, size) focal(x, w = matrix(1/size,
                                               nc=size, nr=size))


  if (parallel) {
    # go parallel
    no_threads <- parallel::detectCores() - free_threads

    ## Initiate cluster
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl, c("sizes", "bin_plant"),
                            environment())
    stack_plant <- parallel::parLapply(cl,
                                       sizes,
                                       function(size) fun(bin_plant, size))

    ## finish
    parallel::stopCluster(cl)

    ## Initiate cluster
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl, c("sizes", "bin_sky"),
                            environment())
    stack_sky <- parallel::parLapply(cl,
                                     sizes,
                                     function(size) fun(bin_sky, size))

    ## finish
    parallel::stopCluster(cl)
  } else {
    stack_plant <- lapply(sizes, function(size) fun(bin_plant, size))
    stack_sky <- lapply(sizes, function(size) fun(bin_sky, size))
  }

  stack_plant <- Map(function(i) stack_plant[[i]] *
                       sizes[i], seq_along(sizes) )
  stack_sky <- Map(function(i) stack_sky[[i]] * sizes[i],
                   seq_along(sizes) )

  # Flatten over sample_size
  fun <- function(x) {
    x[x > sample_size] <- sample_size
    x
  }
  stack_plant <- Map(fun, stack_plant)
  stack_sky <- Map(fun, stack_sky)

  # Create label
  fun <- function(x, size) {x * 100 + max(sizes) - size}
  stack_plant <- Map(fun, stack_plant, sizes)
  stack_sky <- Map(fun, stack_sky, sizes)

  # Find the minimum size that gets at least the sample_size
  flattened_plant <- calc(stack(stack_plant), max)
  ## transform the label into data
  min_size_plant <- flattened_plant - trunc(flattened_plant / 100) * 100
  min_size_plant <- max(sizes) - min_size_plant

  flattened_sky <- calc(stack(stack_sky), max)
  ## transform the label into data
  min_size_sky <- flattened_sky - trunc(flattened_sky / 100) * 100
  min_size_sky <- max(sizes) - min_size_sky
  rm(flattened_plant, flattened_sky, stack_plant,
     stack_sky)

  # Filter out non-mixed-pixels
  mixed_pixel_mask <- !bin_plant & !bin_sky & m
  min_size_plant[!mixed_pixel_mask] <- NA
  min_size_sky[!mixed_pixel_mask] <- NA

  # For unknow reason, I needed to run this.
  min_size_plant <- round(min_size_plant)
  min_size_sky <- round(min_size_sky)

  # Get the required sizes
  sizes_plant <- unique(min_size_plant)
  sizes_sky <- unique(min_size_sky)

  # Remove the sizes that should not be there.
  temp <- match(sizes_plant, sizes)
  for (i in sizes_plant[is.na(temp)]) {
    index <- which.min(abs(sizes - i))
    min_size_plant[min_size_plant == i] <- sizes[index]
  }

  temp <- match(sizes_sky, sizes)
  for (i in sizes_sky[is.na(temp)]) {
    index <- which.min(abs(sizes - i))
    min_size_sky[min_size_sky == i] <- sizes[index]
  }

  # Get the required sizes --again
  sizes_plant <- unique(min_size_plant)
  sizes_sky <- unique(min_size_sky)

  # Filter the blue layer
  blue_plant <- blue_sky <- blue
  blue_plant[!bin_plant] <- NA
  blue_sky[!bin_sky] <- NA

  # Generate the required sizes
  fun <- function(x, y) focal(x, matrix(1, y, y),
                              fun = median, na.rm = TRUE)

  if (parallel) {
    ## Initiate cluster
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl,
                            c("sizes_plant", "blue_plant"),
                            environment())
    stat_plant <- parallel::parLapply(cl,
                                      sizes_plant,
                                      function(y) fun(blue_plant, y))
    ## finish
    parallel::stopCluster(cl)

    ## Initiate cluster
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl,
                            c("sizes_sky", "blue_sky"),
                            environment())
    stat_sky <- parallel::parLapply(cl,
                                    sizes_sky,
                                    function(y) fun(blue_sky, y))

    ## finish
    parallel::stopCluster(cl)
  } else {
    stat_plant <- Map(function(y) fun(blue_plant, y),
                      sizes_plant)
    stat_sky <- Map(function(y) fun(blue_sky, y),
                    sizes_sky)
  }

  # Find the most local stat
  most_local_stat_plant <- most_local_stat_sky <- raster(blue)

  for (i in seq_along(sizes_plant)) {
    indices <- min_size_plant == sizes_plant[i]
    most_local_stat_plant[indices] <- stat_plant[[i]][indices]
  }

  for (i in seq_along(sizes_sky)) {
    indices <- min_size_sky == sizes_sky[i]
    most_local_stat_sky[indices] <- stat_sky[[i]][indices]
  }

  stat_plant <- most_local_stat_plant
  stat_sky <- most_local_stat_sky
  rm(most_local_stat_plant, most_local_stat_sky)

  ## fill the voids
  r <- raster(ncol = 100, nrow = 100)
  extent(r) <- extent(0, ncol(blue), 0, nrow(blue))
  projection(r) <- NA

  fun <- function(stat_x) {
    aux <- resample(stat_x, r)
    aux <- .interpolate_akima(aux)
    i <- 0
    while (i < 3) {
      i <- i + 1
      aux <- focal(aux, matrix(1,3,3), fun = mean,
                   na.rm = TRUE, NAonly = TRUE)
    }
    aux <- resample(aux, stat_x)
    cover(stat_x, aux)
  }

  stat_plant <- fun(stat_plant)
  stat_sky <- fun(stat_sky)

  # Equation 12 from Leblanc et al. (2005)
  # DOI: 10.1016/j.agrformet.2004.09.006
  Leblanc2005_eq12 <- function(DN, DNmin, DNmax) {
    (DN - DNmin) / (DNmax - DNmin)
  }

  # combine pure pixels with mixed pixels
  gf <- Leblanc2005_eq12(blue, stat_plant, stat_sky)
  gf[!mixed_pixel_mask] <- NA
  gf[gf > 0.9] <- 1
  gf[gf < 0] <- 0
  gf[bin_sky] <- 1
  gf[bin_plant] <- 0

  gf
}



#' Subpixel to binary classification
#'
#' Transform a subpixel classification from \code{mblt_subpixel} into a
#' binarized image.
#'
#' This function is useful to export the result of \code{\link{mblt_subpixel}}
#' to a software that does not handle subpixel classification, such as Hemisfer
#' and CIMES.
#'
#' @param subpixel \linkS4class(RasterLayer). The result of a call to
#'   \code{\link{mblt_subpixel}}.
#' @param segmentation \linkS4class(RasterLayer).This should match with the type
#'   of sky grid that will be used in the software next in the processing
#'   pipeline.
#'
#'
subpixel_to_bin <- function (subpixel, segmentation) {
  stopifnot(class(subpixel) == "RasterLayer")
  stopifnot(class(segmentation) == "RasterLayer")
  compareRaster(subpixel, segmentation)

  subpixel[is.na(subpixel)] <- 0

  fun <- function(x) {
    no_of_pixels <- round(mean(x) * length(x))

    if (no_of_pixels > 0 &
        no_of_pixels != length(x)) {
      indices <- order(x, decreasing = TRUE)[1:no_of_pixels]
      x[indices] <- 1
      x[x != 1] <- 0

    } else {
      x <- as.numeric(x > 0.5)
    }
    x
  }

  .cells <- subpixel
  values(.cells) <- 1:ncell(subpixel)
  .cells <- tapply(values(.cells), values(segmentation), function(x) x)
  .cells <- unlist(.cells)
  bin <- tapply(values(subpixel), values(segmentation), fun)
  subpixel[.cells] <- unlist(bin)
  subpixel
}

