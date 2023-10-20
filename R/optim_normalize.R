#' Optimize normalize parameters
#'
#' Wrapper function for [bbmle::mle2()]. Optimize normalize parameters by
#' maximizing [colorfulness()] and minimizing saturation.
#'
#' @inheritParams enhance_caim
#' @inheritParams normalize
#' @inheritParams bbmle::mle2
#'
#' @family Tool Functions
#'
#' @return Numeric vector of length two.
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' m <- !is.na(z)
#' mn_mx <- optim_normalize(caim, m)
#' ncaim <- normalize(caim, mn_mx[1], mn_mx[2], TRUE)
#' plotRGB(ncaim*255)
#' plotRGB(normalize(caim)*255)
#' }
optim_normalize <- function(caim, m, mn = NULL, mx = NULL, method = "BFGS")  {
  terra::compareGeom(caim, m)
  stopifnot(terra::nlyr(caim) == 3)

  names(caim) <- names(read_caim())
  caim[is.na(caim)] <- 0
  if(is.null(mn)) mn <- quantile(caim$Blue[m], 0.01)
  if(is.null(mx)) mx <- quantile(caim$Blue[m], 0.99)

  .get_index <- function(mn, mx) {
    .caim <- normalize(caim, mn, mx, TRUE)
    area <- (sum(.caim$Blue[m] == 0 | .caim$Blue[m] == 1) / sum(m[])) * 100
    cf <- 0
    try(cf <- colorfulness(.caim, m), silent = TRUE)
    (100 - cf) * area
  }
  fit <- bbmle::mle2(.get_index, list(mn = mn, mx = mx), method = method)
  fit@coef
}
