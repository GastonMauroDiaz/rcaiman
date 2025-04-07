#' Optimize a parameter of the function [normalize_minmax()]
#'
#' Wrapper function for [stats::optim()]. Optimize the `mx` argument of the
#' function [normalize_minmax()] by maximizing [colorfulness()] and minimizing
#' saturation.
#'
#' @inheritParams enhance_caim
#' @inheritParams ootb_mblt
#' @inheritParams stats::optim
#'
#' @family Tool Functions
#'
#' @return Numeric vector of length one. The values for using as `mx` argument
#'   with [normalize_minmax()].
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
#' r <- normalize_minmax(caim$Blue, mn, mx, TRUE)
#'
#' bin <- find_sky_pixels(r, z, a)
#' mblt <- ootb_mblt(r, z, a, bin)
#' plot(mblt$bin)
#'
#' mx <- optim_normalize(caim, mblt$bin)
#' ncaim <- normalize_minmax(caim, mx = mx, force_range = TRUE)
#' plotRGB(ncaim*255)
#' plotRGB(normalize_minmax(caim)*255)
#' percentage_of_clipped_highlights(ncaim$Blue, m)
#' }
optim_normalize <- function(caim, bin, method = "BFGS")  {

  terra::compareGeom(caim, bin)
  stopifnot(terra::nlyr(caim) == 3)

  caim[is.na(caim)] <- 0
  mx <- quantile(caim$Blue[bin], 0.99)

  .get_index <- function(mx) {
    .caim <- normalize_minmax(caim, mx = mx, force_range = TRUE)
    f_area <- (sum(.caim$Blue[bin] == 0 | .caim$Blue[bin] == 1) / sum(bin[]))
    tryCatch(cf <- colorfulness(.caim, bin),
             error = function(e) 0)
    (100 - cf) * (f_area * 100)
  }
  opt_result <- stats::optim(mx, .get_index, method = method)
  mx <- opt_result$par %>% unname()
  mx
}
