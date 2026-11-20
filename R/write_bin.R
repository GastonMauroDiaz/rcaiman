#' Write and read binarized images
#'
#' Wrapper functions around [write_caim()] and [read_caim()] to read and write binary masks.
#'
#' @section Functions:
#' \describe{
#'   \item{`write_bin`}{Write a one-layer logical [terra::SpatRaster-class]
#'     to disk as a GeoTIFF (`GTiff`, `INT1U`). No return value.}
#'   \item{`read_bin`}{Read a file with values `255` and/or `0`, such as the one
#'      produced by `write_bin` (see *Details*), and return a logical
#'     [terra::SpatRaster-class] (`TRUE` for `255`, `FALSE` for `0`).}
#' }
#'
#' @param path character vector of length one. File path to read or write. See
#'   examples.
#' @param bin logical [terra::SpatRaster-class] with a single layer.
#'
#' @return See *Functions*
#'
#' @name read_bin
#' @rdname read_bin
#' @aliases write_bin
#'
#' @export
#'
#' @seealso [read_caim()], [write_caim()].
#'
read_bin <- function(path) {
  .check_vector(path, "character", 1)
  .assert_file_exists(path)

  r <- read_caim(path, fill_na = "value", fill_value = 255)
  as.logical(r/255)
}

#' @rdname read_bin
#' @export
write_bin <- function(bin, path) {
  .assert_logical_mask(bin)
  .check_vector(path, "character", 1)

  write_caim(bin * 255, path, 8)
}

