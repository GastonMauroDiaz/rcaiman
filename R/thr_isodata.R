#' Threshold calculated with the isodata method
#'
#' Threshold calculated with the algorithm by
#' \insertCite{isodata;textual}{rcaiman}, which was recommended by
#' \insertCite{Jonckheere2005;textual}{rcaiman}.
#'
#' The implementation is based on
#' \href{http://fiji.sc/Auto_Threshold#IsoData}{the IsoData method of Auto
#' Threshold ImageJ plugin by Gabriel Landini}, which is now available in the
#' 'autothresholdr' package (\code{\link[autothresholdr]{auto_thresh}}).
#' However, I found this implementarion more versatile since it is not
#' restricted to an 8-bit input.
#'
#' @param x Numeric vector or a single-column \emph{matrix} or \emph{data.frame}
#'   able to be coerced to numeric.
#'
#' @return Numeric vector of length one.
#' @export
#'
#' @family Binarization Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' caim <- read_caim()
#' r <- gbc(caim$Blue)
#' thr <- thr_isodata(values(r))
#' bin <- apply_thr(r, thr)
#' plot(bin)
thr_isodata <- function(x) {
  if(is(x, "matrix") | is(x, "data.frame")) {
    x <- as.numeric(x)
  }
  stopifnot(is(x, "numeric"))
  x <- x[!is.na(x)]
  size <- length(x)
  if (size <= 1) stop("length(x) must be greater than 1.")
  if (size > 500) size <- 500
  if (stats::sd(sample(x, size)) == 0) {
    thr <- x[1]
    warning("sd(x) should be greater than 0.")
  } else {
    thr <- mean(x)
    thr.back <- 0
    while (thr != thr.back) {
      thr.back <- thr
      x0 <- x[x <= thr]
      x1 <- x[x > thr]
      thr <- (mean(x0) + mean(x1)) / 2
    }
  }
  thr
}
