#' Optimize normalize parameters
#'
#' Wrapper function for [bbmle::mle2()]. Optimize normalize parameters by
#' maximizing [colorfulness()] and minimizing saturation.
#'
#' @inheritParams enhance_caim
#' @inheritParams ootb_mblt
#' @inheritParams bbmle::mle2
#'
#' @family Tool Functions
#'
#' @return Numeric vector of length one. The values for using as `mx` argument
#'   with [normalize()].
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#'
#' mn <- quantile(caim$Blue[m], 0.01)
#' mx <- quantile(caim$Blue[m], 0.99)
#' r <- normalize(caim$Blue, mn, mx, TRUE)
#'
#' bin <- find_sky_pixels(r, z, a)
#' mblt <- ootb_mblt(r, z, a, bin)
#' plot(mblt$bin)
#'
#' mx <- optim_normalize(caim, mblt$bin)
#' ncaim <- normalize(caim, mx = mx, force_range = TRUE)
#' plotRGB(ncaim*255)
#' plotRGB(normalize(caim)*255)
#' percentage_of_clipped_highlights(ncaim$Blue, m)
#' }
optim_normalize <- function(caim, bin, method = "BFGS")  {
  terra::compareGeom(caim, bin)
  stopifnot(terra::nlyr(caim) == 3)

  names(caim) <- names(read_caim())
  caim[is.na(caim)] <- 0
  mx <- quantile(caim$Blue[bin], 0.99)

  .get_index <- function(mx) {
    .caim <- normalize(caim, mx = mx, force_range = TRUE)
    area <- (sum(.caim$Blue[bin] == 0 | .caim$Blue[bin] == 1) / sum(bin[])) * 100
    cf <- 0
    try(cf <- colorfulness(.caim, bin), silent = TRUE)
    (100 - cf) * area
  }
  fit <- bbmle::mle2(.get_index, list(mx = mx), method = method)
  fit@coef
}
