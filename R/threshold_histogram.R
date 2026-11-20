# Gaussian smoothing of equally spaced histogram frequencies
.smooth_hist_gaussian <- function(h, sigma, pad = c("reflect","replicate")) {
  if (is.null(sigma)) {
    ys <- h
  } else {
    pad <- match.arg(pad)
    stopifnot(length(h) >= 2, is.finite(sigma), sigma > 0)
    r <- ceiling(3 * sigma)
    x <- -r:r
    k <- exp(-(x*x) / (2 * sigma * sigma))
    k <- k / sum(k)  # normalize kernel to unit area

    # build padded vector
    if (pad == "replicate") {
      padL <- rep(h[1], r); padR <- rep(h[length(h)], r)
    } else { # reflect
      reflect <- function(i) {
        period <- 2 * (length(h) - 1)
        j <- (i - 1) %% period
        h[1 + pmin(j, period - j)]
      }
      padL <- reflect(seq.int(1 - r, 0))
      padR <- reflect(length(h) + seq_len(r))
    }
    hp <- c(padL, h, padR)

    # convolution (centered, same length)
    ys <- as.numeric(stats::filter(hp, k, sides = 2))
    ys <- ys[(r + 1):(length(ys) - r)]
  }
  ys
}


.rosin_corner <- function(h, peak_dn, slope_reduction = FALSE) {
  # Coordinates are bin positions, starting at zero; the endpoint is one
  # bin beyond the histogram. No observation rounding takes place here.
  xs <- seq_along(h) - 1L
  x1 <- peak_dn
  y1 <- h[peak_dn + 1L]
  if (slope_reduction && y1 > mean(h)) y1 <- mean(h)
  x2 <- length(h)
  dx <- x2 - x1
  dy <- -y1
  line <- y1 * (x2 - xs) / dx
  valid <- which(xs > x1 & h < line)
  if (!length(valid)) {
    stop("No corner found beyond the peak on the selected histogram side.", call. = FALSE)
  }
  distances <- (line - h) * dx / sqrt(dx^2 + dy^2)
  i <- valid[which.max(distances[valid])]
  x0 <- xs[i]
  y0 <- h[i]
  t <- ((x0 - x1) * dx + (y0 - y1) * dy) / (dx^2 + dy^2)
  list(x = xs, y = h, p1 = c(x1, y1), p2 = c(x2, 0),
       p0 = c(x0, y0), perpendicular_foot = c(x1 + t * dx, y1 + t * dy),
       max_dist = distances[i])
}

.histogram_corner <- function(h, peak, tail, slope_reduction = FALSE) {
  n <- length(h$counts)
  counts <- h$counts / max(h$counts) * (n - 1L)
  if (tail == "left") {
    counts <- rev(counts)
    peak <- n + 1L - peak
  }
  corner <- .rosin_corner(counts, peak - 1L, slope_reduction)
  corner$index <- as.integer(if (tail == "left") n - corner$p0[1] else corner$p0[1] + 1L)
  corner
}

.plot_histogram_corners <- function(h, corners, tails, title) {
  n <- length(h$counts)
  graphics::plot(seq_len(n) - 1L, h$counts / max(h$counts) * (n - 1L),
    type = "l", asp = 1, xlim = c(-1, n), ylim = c(0, n - 1L), axes = FALSE,
    xlab = "Value", ylab = "Frequency relative to peak", main = title)
  ticks <- seq(0, n - 1L, length.out = 4)
  labels <- h$mids[1] + ticks * diff(h$mids)[1]
  graphics::axis(1, at = ticks, labels = signif(labels, 4))
  frequency_ticks <- seq(0, 1, by = .25)
  graphics::axis(2, at = frequency_ticks * (n - 1L), labels = frequency_ticks)
  for (i in seq_along(corners)) {
    corner <- corners[[i]]
    display_x <- function(v) if (tails[i] == "left") n - 1L - v else v
    col <- c("black", "blue")[(i - 1L) %% 2L + 1L]
    graphics::segments(display_x(corner$p1[1]), corner$p1[2],
                       display_x(corner$p2[1]), corner$p2[2], col = col)
    graphics::points(display_x(corner$p0[1]), corner$p0[2], pch = 19, col = col)
    graphics::segments(display_x(corner$p0[1]), corner$p0[2],
      display_x(corner$perpendicular_foot[1]), corner$perpendicular_foot[2], lty = 2, col = col)
    graphics::abline(v = display_x(corner$p0[1]), lty = 3, col = col)
  }
}
