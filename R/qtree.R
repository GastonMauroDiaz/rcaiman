#' Quad-tree segmentation
#'
#' The quad-tree segmentation algorithm is a top-down process that makes
#' recursive divisions in four equal parts until a condition is satisfied and
#' stops locally. This is the usual implementation of the quad-tree algorithm,
#' so it produces squared segments of different sizes. This particular
#' implementation allows up to five sizes.
#'
#' The algorithm starts splitting the entire image into large squared segments
#' following, depending on the aspect ratio, grids going from \eqn{4 \times 4}
#' to \eqn{1 \times 4} or \eqn{4 \times 1}; then, splits each segment into four
#' sub-segments and calculates the standard deviation of the pixels from
#' \code{r} delimited by each of those segments. The splitting process stops
#' locally if the sum of the standard deviation of the sub-segments minus the
#' standard deviation of the parent segment (named \emph{delta}) is less or
#' equal than the \code{scale_parameter}. If \code{r} has more than one layer,
#' \emph{delta} is calculated separately and \emph{delta} mean is used to
#' evaluate the stopping condition.
#'
#' @inheritParams polar_qtree
#'
#' @return A single layer image of the class \linkS4class{SpatRaster} with
#'   integer values.
#'
#' @export
#'
#' @family Segmentation Functions
#'
#' @examples
#' \donttest{
#' caim <- read_caim()
#' plot(caim)
#' caim <- normalize(caim, 0, 255)
#' seg <- qtree(caim, scale_parameter = 0.5)
#' plot(caim$Blue)
#' plot(extract_feature(caim$Blue, seg))
#' plot(extract_feature(seg, seg, length))
#' }
qtree <- function(r, scale_parameter = 0.2) {
  stopifnot(class(r) == "SpatRaster")
  stopifnot(is.numeric(scale_parameter))
  stopifnot(length(scale_parameter) == 1)
  r <- terra::deepcopy(r)

  x <- round(min(ncol(r), nrow(r))/4)
  y <- paste0(c(1, rep(0, (nchar(x)-1))), collapse = "") %>% as.numeric()
  x <- (substr(x,1,1) %>% as.numeric())
  if (!.is_even(x)) x <- x + 1
  x <- x*y
  wds <- x / 2^(1:10)
  wds <- wds[wds > 2^2]
  wds <- wds[.is_whole(wds)]
  if (length(wds) > 5) wds <- wds[1:5]
  if (length(wds) < 1) stop("A greater number of pixels is required")

  e <- terra::ext(0, ((trunc(ncol(r)/wds[1])) + 1) * wds[1],
                  0, ((trunc(nrow(r)/wds[1])) + 1) * wds[1])
  .r <- terra::extend(r, e)

  ges <- Map(chessboard, .r, wds)
  names(ges) <- NULL

  .calc_delta_single_layer <- function(r) {
    .sd <- function(x) {
      x <- sd(x, na.rm = TRUE)
      if (is.na(x)) x <- 0
      x
    }
    sd_now <- Map(function(g) extract_feature(r, g, .sd,
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
    delta <- Map(.calc_delta_single_layer, as.list(.r))
    delta <- Map(function(i) {
      x <- Map(function(j) delta[[j]][[i]], seq_along(delta))
      apply(as.data.frame(x), 1, mean)
    }, seq_along(delta[[1]]))
  } else {
    delta <- .calc_delta_single_layer(.r)
  }

  .it_should_be_splited <- function(delta) {
    it_should_be_splited <- Map(function(x) x > scale_parameter, delta)
    it_should_be_splited <- Map(function(i) {
      terra::subst(ges[[i]],
                   names(delta[[i]]) %>%
                     as.numeric(),
                   it_should_be_splited[[i]])
    }, seq_along(it_should_be_splited))
    ges <- ges[-1]
    cte <- rep(0, nchar(max(terra::rast(ges)[]))) %>% c(1, .) %>%
      paste0(., collapse = "") %>% as.numeric()
    ges <- Map(function(i) {
      r <- ges[[i]]
      r[r == 0] <- NA
      r + i*cte
    }, seq_along(ges))
    ges <- terra::rast(ges)
    it_should_be_splited <- terra::rast(it_should_be_splited)
    seg <- ges * it_should_be_splited
    seg <- max(seg)
    size <- extract_feature(seg, seg, length)
    i <- unique(size[])
    i <- i[!(sqrt(i) %in% wds)]
    i <- i[!is.na(i)]
    if (length(i) > 0) {
      for (u in i) {
        seg[size == u] <- 0
      }
    }
    seg
  }
  seg <- .it_should_be_splited(delta)
  m <- seg == 0
  foo <- terra::rast(ges) * m
  # cte <- paste0(c(1, rep(0, nchar(max(foo[])))), collapse = "") %>% as.numeric()
  foo <- Map(function(i) {
               m <- extract_feature(foo[[i]], foo[[i]], length) == wds[[i]]^2
               m * (length(wds) - i)
             },
             seq_along(wds))
  foo <- length(ges) - max(terra::rast(foo), na.rm = TRUE)
  foo[is.na(foo)] <- 0
  for (i in seq_along(ges)){
    indices <- foo == i
    if (any(indices[] %>% as.logical())) {
      if (i == 1) {
        foo[indices] <- ges[[1]][indices] + length(ges)
      } else {
        foo[indices] <- ges[[i]][indices] + fact
      }
      fact <- max(ges[[i]][indices])
    } else {
      if (i == 1) fact <- 1
    }
  }
  foo <- foo + max(seg[])
  seg <- c(foo*m, seg)
  seg <- max(seg)

  seg <- terra::crop(seg, r)
  names(seg) <- "Quad-tree"
  seg
}
