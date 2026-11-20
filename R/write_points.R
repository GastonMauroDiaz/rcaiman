#' Read and write image points
#'
#' @description
#' Read and write continuous image coordinates stored as CSV files.
#'
#' @section Functions:
#' \describe{
#'   \item{`read_points`}{Read a CSV file and return a `data.frame` with columns
#'     `"x"` and `"y_down"`. Additional columns are preserved.}
#'   \item{`write_points`}{Write continuous image coordinates to disk as a CSV
#'     file. No return value.}
#' }
#'
#' @param path character vector of length one. File path to read or write.
#' @param sampling_coord `data.frame` with columns `"x"` and `"y_down"`,
#'   representing continuous image coordinates. Additional columns are
#'   preserved.
#'
#' @details
#' Coordinates use an image reference system with its origin at the upper-left
#' corner. The `"x"` coordinate increases to the right and `"y_down"` increases
#' downward. The center of the upper-left cell is `(0.5, 0.5)`.
#'
#' `write_points()` writes CSV files. If `path` has no extension, `".csv"` is
#' added. It places `"x"` and `"y_down"` first, in that order. This CSV schema
#' defines the point-data interoperability contract between \pkg{rcaiman} and
#' [rcaiman-projects](https://github.com/GastonMauroDiaz/rcaiman-projects).
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
#' @seealso [indices_from_coord()] to assign coordinates to raster cells and
#'   [coord_from_indices()] to obtain cell-center coordinates.
#'
#' @examples
#' \dontrun{
#' sampling_coord <- data.frame(
#'   x = c(10.2, 20.8),
#'   y_down = c(15.7, 25.3),
#'   point_id = c("p1", "p2")
#' )
#' my_file <- file.path(tempdir(), "pts.csv")
#' write_points(sampling_coord, my_file)
#' sampling_coord_from_disk <- read_points(my_file)
#' sampling_coord_from_disk
#' }
read_points <- function(path) {
  .check_vector(path, "character", 1)
  .assert_file_exists(path)
  sampling_coord <- utils::read.csv(path)
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down")
  )
  sampling_coord
}

#' @rdname read_points
#' @export
write_points <- function(sampling_coord, path) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down")
  )
  .check_vector(path, "character", 1)

  file_name <- basename(path)
  file_name <-  .assert_file_extension(file_name, "csv")

  nms <- names(sampling_coord)
  other_columns <- setdiff(nms, c("x", "y_down"))
  sampling_coord <- sampling_coord[, c("x", "y_down", other_columns)]

  utils::write.csv(
    sampling_coord,
    file.path(dirname(path), file_name),
    row.names = FALSE
  )
}

