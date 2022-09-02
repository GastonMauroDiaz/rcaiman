qtree <- function(r, scale_parameter = 0.2) {
  stopifnot(class(r) == "SpatRaster")

  x <- round(min(ncol(r), nrow(r))/4)
  y <- paste0(c(1, rep(0, (nchar(x)-1))), collapse = "") %>% as.numeric()
  x <- (substr(x,1,1) %>% as.numeric())
  if (!.is_even(x)) x <- x + 1
  x <- x*y
  wds <- x / 2^(1:10)
  wds <- wds[wds > 2^2]
  wds <- wds[.is_whole(wds)]
  if (length(wds) < 1) stop("A greater number of pixels is required")

  .r <- terra::subset(r, 1)
  .ges <- Map(chessboard, .r, wds)
  m <- Map(function(i) extract_feature(.r, .ges[[i]], length) == wds[[i]]^2,
           seq_along(wds))
  rm(.r)
  ges <- Map(function(i) .ges[[i]] * m[[i]], seq_along(wds))

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
    delta <- Map(.calc_delta_single_layer, as.list(r)) %>% suppressWarnings()
    delta <- Map(function(i) {
      x <- Map(function(j) delta[[j]][[i]], seq_along(delta))
      apply(as.data.frame(x), 1, mean)
    }, seq_along(delta[[1]]))
  } else {
    delta <- .calc_delta_single_layer(r) %>% suppressWarnings()
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
  seg <- max(seg)
  seg[seg == 0] <- NA
  middle <- length(wds) / 2
  middle <- round(middle)
  seg <- terra::cover(seg, .ges[[middle]])
  names(seg) <- "Quad-tree"
  seg
}
