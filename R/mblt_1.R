#' Threshold function
#'
#' Threshold value in function of the background digital number.
#'
#' @param dn Numeric vector. Digital number of the background. The background
#'   should be lighter than the objects. In canopy photography, the typical
#'   background is the sky.
#' @param w Numeric vector of length one. Weight. See reference.
#' @param type Character vector of length one. Default is "Generic". Currently,
#'   the only available calibrated values are from
#'   \insertCite{Diaz2018}{rcaiman}. Use "Nikon_Coolpix_5700" to use them.
#' @param intercept,slope Numeric vector of length one. Default is NULL. You can
#'   provide your own calibrated values. See \insertCite{Diaz2018}{rcaiman} for
#'   details.
#'
#' @export
#'
#' @references
#' \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' thr_fun(125)
thr_fun <- function (dn,
                     w = 0.5,
                     type = "Generic",
                     intercept = NULL,
                     slope = NULL) {
  stopifnot(class(dn) == "numeric")
  stopifnot(length(w) == 1)
  stopifnot(class(w) == "numeric")
  stopifnot(class(type) == "character")
  if (!is.null(intercept)) stopifnot(length(intercept) == 1)
  if (!is.null(slope)) stopifnot(length(slope) == 1)

  if (all(is.null(intercept), is.null(slope))) {

    my_coef <- switch(type,
                      Generic = c(-8, 1),
                      Nikon_Coolpix_5700 = c(-7.7876, 0.9485))

    intercept <- my_coef[1]
    slope <- my_coef[2]
  }

  if (length(dn[dn > 255]) > 0) {dn[dn > 255] <- 255}
  thr <- intercept  + slope * dn * w
  if (length(thr[thr < 0]) > 0) {thr[thr < 0] <- 0}
  thr
}







