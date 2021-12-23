#' @import methods
#' @import raster
#' @importFrom magrittr %>%
#' @importFrom stats lm poly coefficients median IQR sd splinefun
#' @importFrom colorspace sRGB
NULL

# https://github.com/tidyverse/magrittr/issues/29
#' @importFrom utils globalVariables
NULL
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
