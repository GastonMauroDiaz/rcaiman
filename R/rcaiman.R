#' rcaiman: An R package for CAnopy IMage ANalysis
#'
#' Solutions for processing and pre-processing canopy images, particularly
#' hemispherical photographs, including non-circular ones, such as certain
#' pictures taken with auxiliary fisheye lens attached to smartphones.
#'
#' @section Flagship functions:
#' \itemize{
#'
#' \item \code{\link{expand_noncircular}}
#' \item \code{\link{fisheye_to_equidistant}}
#' \item \code{\link{ootb_mblt}}
#' \item \code{\link{ootb_obia}}
#' }
#'
#' @section Batch Processing:
#'
#'   Batch processing can be easily performed with standard R programming. Below
#'   is an example that can be used as a template.
#'
#'   \preformatted{
#' require(rcaiman)
#'
#' input_folder <- "c:/Users/janedoe/pics/"
#' output_folder <- "c:/Users/janedoe/bins/"
#' files <- dir(input_folder, full.names = TRUE)
#'
#' for (i in 1:length(files)) {
#'  caim <- read_caim(file.path(files[i]))
#'  blue <- gbc(caim$Blue)
#'  bin <- apply_thr(blue, thr_isodata(blue[]))
#'  write_bin(bin, file.path(output_folder, basename(files[i])))
#' }
#'   }
#'
#' @docType package
#' @name rcaiman
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
