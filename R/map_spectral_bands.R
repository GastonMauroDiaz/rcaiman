#' Map spectral bands
#'
#' @description
#' Derive named spectral bands from CFA or sensor bands using algebraic
#' expressions.
#'
#' @details
#' The input `caim` is typically obtained from [split_cfa()].
#'
#' `spectral_mapping` defines how target bands are computed from CFA channels
#' using algebraic expressions. Each component must be a character vector of
#' length one representing an expression involving:
#'
#' - CFA channel names present in `caim`
#' - numeric constants
#' - the operators `+`, `-`, `*`, `/`, `(`, and `)`
#'
#' The names of `spectral_mapping` determine the names of the output layers.
#'
#' Example:
#'
#' ```
#' spectral_mapping <- list(
#'   Red = "(Magenta + Yellow)/2",
#'   Green = "Green",
#'   Blue = "Cyan"
#' )
#' ```
#'
#' Expressions are evaluated in an environment where channel names correspond
#' to layer names of `caim`.
#'
#' @param caim numeric [terra::SpatRaster-class]. Layers must have unique names.
#' @param spectral_mapping named list of character vectors defining algebraic
#' expressions used to derive target spectral bands from CFA channels.
#'
#' @return
#' Numeric [terra::SpatRaster-class] with one layer per element of
#' `spectral_mapping`. Layer names match `names(spectral_mapping)`.
#' The output preserves the geometry of `caim`.
#'
#' @export
#' @examples
#' caim <- read_caim()
#' map_spectral_bands(caim, list(Red = "(Magenta + Yellow)/2",
#'                            Green = "Green",
#'                            Blue = "Cyan"))
map_spectral_bands <- function(caim, spectral_mapping) {
  .assert_spatraster(caim)
  .check_spectral_mapping(spectral_mapping, names(caim))

  if (length(names(caim)) != length(unique(names(caim))))
    stop("Layer names of `caim` must be unique.")


  bands <- as.list(caim)
  names(bands) <- names(caim)

  env <- list2env(bands, parent = baseenv())

  out <- vector("list", length(spectral_mapping))
  for (i in seq_along(spectral_mapping)) {

    expr <- parse(text = spectral_mapping[[i]])[[1]]

    out[[i]] <- eval(expr, env)

  }

  out <- terra::rast(out)
  names(out) <- names(spectral_mapping)
  out
}
