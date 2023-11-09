#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' @import methods
#' @import terra
#' @import filenamer
#' @importFrom magrittr %>%
#' @importFrom stats lm poly coefficients IQR sd splinefun
#' @importFrom colorspace sRGB
#' @importClassesFrom lidR LAS
NULL

# https://groups.google.com/g/rdevtools/c/qT6cJt6DLJ0
# spurious importFrom to avoid note
#' @importFrom Rdpack c_Rd
NULL


# https://github.com/tidyverse/magrittr/issues/29
#' @importFrom utils globalVariables
NULL
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
