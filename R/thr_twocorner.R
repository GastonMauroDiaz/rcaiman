#' Calculate the thresholds of the two-corner method
#'
#' @inheritParams thr_isodata
#' @param diagnose Logical vector of length one. If `TRUE`, then a figure will
#'   be send to the graphical device showing the geometrical contruction made to
#'   find the threholds.
#' @param window_length Numeric vector. One or many odd numbers indicating the
#'   size of the window to be used to smooth the histogram.
#' @param slope_reduction Logical vector of length one. If `TRUE`, then the
#'   slope reduction method by Macfarlane (2011) will be applied.
#'
#' @returns a list.
#' @export
#'
#' @examples
#' \dontrun{
#' }
thr_twocorner <- function(x,
                          diagnose = FALSE,
                          window_length = NULL,
                          slope_reduction = TRUE) {

  .find_dn_of_first_empty_bin <- function(h, peak) {
    # “the line starts at the largest bin and finishes at the first empty bin of
    # the histogram following the last filled bin.” ([Rosin, 2001, p. 2084]
    h[1:peak] <- max(h)
    zero_bins <- h == 0
    if (any(zero_bins)) {
      x2 <- which.max(zero_bins) - 1
    } else {
      x2 <- 256
    }
    x2
  }

  .moving_average <- function(y, window_length = 5) {
    if (!is.null(window_length)) {
      if (window_length %% 2 == 0) stop("Window length must be odd")
      pad <- floor(window_length / 2)
      y_padded <- c(rep(y[1], pad), y, rep(y[length(y)], pad))
      y <- stats::filter(y_padded,
                         rep(1 / window_length, window_length), sides = 2)
      y <- y[!is.na(y)]
    }
    y
  }

  .exact_geometric_construction <- function(h, peak_dn, window_length = NULL,
                                            slope_reduction = FALSE) {
    # Rosin, P. L. (2001). Unimodal thresholding. Pattern Recognition, 34(11),
    # 2083–2096. https://doi.org/10.1016/S0031-3203(00)00136-9

    # p2
    x2 <- .find_dn_of_first_empty_bin(h, peak_dn)
    y2 <- 0

    h <- .moving_average(h, window_length)

    if (slope_reduction && h[peak_dn + 1] > mean(h)) {
      # Macfarlane, C. (2011). Classification method of mixed pixels does not
      # affect canopy metrics from digital images of forest overstorey.
      # Agricultural and Forest Meteorology, 151(7), 833–840.
      # https://doi.org/10.1016/j.agrformet.2011.01.019

      # p1
      x1 <- peak_dn
      y1 <- mean(h)
    } else {
      # p1
      x1 <- peak_dn
      y1 <- h[peak_dn + 1]
    }

    # Prepara las coordenadas de todos los puntos del histograma
    # (xi, yi) con xi en 0:255
    xs <- 0:255
    ys <- h

    # Parámetros de la línea p1→p2
    dx <- x2 - x1
    dy <- y2 - y1
    denom <- sqrt(dx^2 + dy^2)

    # Función de distancia perpendicular punto→línea
    perp_dist <- function(x0, y0, x1, y1, dx, dy, denom) {
      num <- abs(dy * x0 - dx * y0 + x2 * y1 - y2 * x1)
      num / denom
    }

    # Vectorizamos
    dists <- perp_dist(
      x0 = xs,
      y0 = ys,
      x1 = x1,
      y1 = y1,
      dx = dx,
      dy = dy,
      denom = denom
    )

    # perpendiculares incluidas entre p1 y p2 y sobre la linea p1→p2
    l <- line(c(x1,x2), c(y1, y2))
    if (any(is.na(coef(l)))) stop("Peak detection failed")
    lfn <- function(x) coef(l)[1] + coef(l)[2] * x
    valid <- (xs + 1)[ys < lfn(xs) & xs > x1]

    if (length(valid) == 0) stop("Empty valid range")

    # punto con máxima distancia entre los válidos
    i_max <- valid[which.max(dists[valid])]
    max_dist <- dists[i_max]

    # Proyección escalar t = ((x0-x1)*dx + (y0-y1)*dy) / (dx^2+dy^2)
    # Parametro t de la proyección (p1 + t*(p2-p1) es el pie)
    x0 <- xs[i_max]
    y0 <- ys[i_max]
    t <- ((x0 - x1)*dx + (y0 - y1)*dy) / (dx^2 + dy^2)

    # Coordenadas del punto proyectado (pie de la perpendicular)
    xp <- x1 + t * dx
    yp <- y1 + t * dy

    # corner
    list(
      x = xs,
      y = ys,
      p1 = c(x1, y1),
      p2 = c(x2, y2),
      p0 = c(x0, y0),
      perpendicular_foot = c(xp, yp),
      max_dist = max_dist
    )
  }

  .revert_normalization <- function(x, mn, mx) x * (mx - mn) + mn

  .fun <- function(window_length) {
    mn <- min(x)
    mx <- max(x)
    vals <- (normalize_minmax(x, mn, mx) * 255) %>% round()
    h <- as.numeric(table(factor(vals, levels = 0:255)))

    # To draw in a square
    h <- h/max(.moving_average(h, window_length))*255

    # 2) Find peaks (Macfarlane 2011)
    DNL1 <- 5; DNL2 <- 55; DNR1 <- 200; DNR2 <- 250

    repeat {
      DNMAX_left  <- which.max(h[DNL1:DNL2]) + (DNL1 - 1)
      ok_left  <- (DNL2 - DNMAX_left  >= 10)
      if (ok_left) break
      if (!ok_left && DNL2 < 255) DNL2 <- min(255, DNL2 + 25)
      if (DNL2 >= 255) {
        stop("Left peack close to the maximum bin")
      }
    }

    repeat {
      DNMAX_right <- which.max(h[DNR1:DNR2]) + (DNR1 - 1)
      # Condiciones de aceptación para cada pico
      ok_right <- (DNMAX_right - DNR1 >= 10)
      if (ok_right) break
      if (!ok_right && DNR1 > 0) DNR1 <- max(0,   DNR1 - 25)
      if (DNR1 <= 0) {
        stop("Right peak close to the minimum bin")
      }
    }

    if ((DNMAX_right - DNMAX_left) < 20) stop("Unimodal histogram")


    # 3) Rosin (2001) modified by Macfarlane (2011)
    left <- .exact_geometric_construction(h, DNMAX_left,
                                          window_length, slope_reduction)
    right <- .exact_geometric_construction(rev(h), 255 - DNMAX_right,
                                           window_length, slope_reduction)

    if (diagnose) {
      .draw <- function(col = "black") {
        segments(x1, y1, x2, y2, lwd=1, col = col)
        points(x1, y1, pch=20, col = col)
        points(x2, y2, pch=20, col = col)
        points(x0, y0, pch=20, col = col)
        points(xp, yp,  pch=20, col = col)
        segments(x0, y0, xp, yp, lty=2, col = col)
        abline(v = x0, lty = 3, col = col)
      }
      xs <- left$x
      ys <- left$y
      x0 <- left$p0[1]
      y0 <- left$p0[2]
      xp <- left$perpendicular_foot[1]
      yp <- left$perpendicular_foot[2]
      x1 <- left$p1[1]
      y1 <- left$p1[2]
      x2 <- left$p2[1]
      y2 <- left$p2[2]

      par_old <- par()
      par(mar = c(0,0,0,0))
      xlim <- c(0,255)
      plot(xs, ys, type='l', asp = 1, xlim = xlim, ylim = xlim,
           axes = FALSE, xlab = "", ylab = "")
      segments(0,-5,255,-5)
      par <- par_old

      .draw()

      xs <- 255- right$x
      ys <- right$y
      x0 <- 255 - right$p0[1]
      y0 <- right$p0[2]
      xp <- 255 - right$perpendicular_foot[1]
      yp <- right$perpendicular_foot[2]
      x2 <- 255 - right$p1[1]
      y2 <- right$p1[2]
      x1 <- 255 - right$p2[1]
      y1 <- right$p2[2]

      .draw("blue")
    }

    lc <- .revert_normalization(left$p0[1]/255, mn, mx)
    uc <- .revert_normalization((255 - right$p0[1])/255, mn, mx)

    list(lc = lc,
         uc = uc,
         tl = lc + (uc - lc) * 0.25,
         tm = lc + (uc - lc) * 0.5,
         th = lc + (uc - lc) * 0.75)

  }

  if (length(window_length) > 1) {
    i <- 1
    repeat {
      if (is.na(window_length[i])) break
      thr <- tryCatch(.fun(window_length[i]), error = function(e) NULL)
      if (!is.null(thr)) break
      i <- i + 1
    }
    if (is.null(thr)) stop("The two-corner method fails after many tries.")
  } else {
    thr <- .fun(window_length)
  }
  thr
}
