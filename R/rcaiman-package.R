#' @section The rcaiman suite:
#' This package is the core toolbox for preprocessing and analysis of canopy
#' photographs. It can be used independently of the other suite components:
#'
#' - [rcaiman.registry](https://github.com/GastonMauroDiaz/rcaiman.registry)
#'   provides calibration helpers and manages equipment specifications.
#' - [rcaiman.registry.data](https://github.com/GastonMauroDiaz/rcaiman.registry.data)
#'   distributes calibration records for reuse.
#' - [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects)
#'   provides a graphical workbench for projects, manual inputs, provenance,
#'   and reusable workflow templates.
#'
#' @section Interface stability from v3:
#' Version 3 introduces a deliberate break with earlier interfaces. From the
#' v3 release onward, established public interfaces remain compatible within
#' each major version. Replacements will be announced with an alternative and
#' a deprecation period before removal in a later major release.
#'
#' Stability is a commitment to interface compatibility; implementations can
#' continue to improve.
#'
#' Functions and model options marked experimental may change without a
#' deprecation period.
#' Lifecycle labels describe interface stability; scientific assumptions and
#' limitations are explained in each function's help. Changes are recorded in
#' the package NEWS.
#'
#' The documentation supports researchers working directly in R or with
#' assistance from language models. Use the documentation for the installed
#' version when preparing scripts.
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
#' @importFrom Rcpp sourceCpp
#' @useDynLib rcaiman, .registration = TRUE
## usethis namespace: end
NULL

#' @import methods
#' @import terra
#' @importFrom magrittr %>%
#' @importFrom magrittr subtract
#' @importFrom stats lm poly coefficients IQR sd splinefun
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom graphics abline par segments
#' @importFrom stats line
NULL

# https://groups.google.com/g/rdevtools/c/qT6cJt6DLJ0
# spurious importFrom to avoid note
#' @importFrom Rdpack c_Rd
NULL


# https://github.com/tidyverse/magrittr/issues/29
#' @importFrom utils globalVariables
NULL
utils::globalVariables(c("."))
