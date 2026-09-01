#' Write and read image points
#'
#' @description
#' Read and write image point coordinates stored as CSV files.
#'
#' @section Functions:
#' \describe{
#'   \item{`write_points`}{Write a `data.frame` with image point coordinates to
#'     disk as a CSV file. No return value.}
#'   \item{`read_points`}{Read a CSV file and return a `data.frame` with, but not exclusively,
#'     columns `"row"`, `"col"`.}
#' }
#'
#' @param path character vector of length one. File path to read or write.
#' @param pts `data.frame` with at least columns `"row"` and `"col"`.
#'
#' @details
#' Point coordinates are stored in raster coordinates, with `"row"` and `"col"`
#' referring to image row and column positions. Additional columns are preserved.
#'
#' `write_points()` writes CSV files. If `path` has no extension, `".csv"` is
#' added. Files written by this function follow the standard rcaiman conventions
#' and can be open by [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects).
#' These conventions include reordering columns when needed, so `"row"` and `"col"`
#' are the first and second columns, respectively.
#'
#'
#' @note
#' File extensions are recognized using [tools::file_ext()].
#'
#' @return See *Functions*
#'
#' @name read_points
#' @rdname read_points
#' @aliases write_points
#'
#' @export
#'
#' @seealso [fibonacci_points()]
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(1000, equidistant_coef)
#' a <- azimuth_image(z)
#' pts <- fibonacci_points(z, a, 10)
#' head(pts)
#' my_file <- file.path(tempdir(), "pts.csv")
#' write_points(pts, my_file)
#' pts_from_disk <- read_points(my_file)
#' head(pts_from_disk)
#' }
read_points <- function(path) {
  .check_vector(path, "character", 1)
  .assert_file_exists(path)
  pts <- read.csv(path)
  .check_columns(pts, required_columns = c("row", "col"))
  pts
}

#' @rdname read_points
#' @export
write_points <- function(pts, path) {
  .check_columns(pts, required_columns = c("row", "col"))
  .check_vector(path, "character", 1)

  file_name <- basename(path)
  file_name <-  .assert_file_extension(file_name, "csv")

  nms <- names(pts)
  if (length(nms) == 2) {
    pts <- pts[, c("row", "col")]
  } else {
    i_row <- grep("row",nms)
    i_col <- grep("col", nms)
    pts <- pts[, c("row", "col", nms[-c(i_row, i_col)])]
  }

  write.csv(pts, file.path(dirname(path), file_name), row.names = FALSE)
}

