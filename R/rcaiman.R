#' rcaiman: An R package for CAnopy IMage ANalysis
#'
#' The rcaiman package provides three categories of important functions:
#' base, diffuse light, and non-difuse. (TODO: revise this categories)
#'
#' Always err on the side of caution, and simplicity. It’s easier to give people
#' more functionality than it is to take away stuff they’re used to
#' https://r-pkgs.org/namespace.html#namespace-workflow
#'
#' @section Base functions:
#' The foo functions ...
#'
#' @docType package
#' @name rcaiman
NULL

#' @import methods
#' @import raster
#' @importFrom magrittr %>%
#' @importFrom stats lm poly coefficients
NULL



# https://github.com/tidyverse/magrittr/issues/29
#' @importFrom utils globalVariables
NULL
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
