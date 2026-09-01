#' Title
#'
#' @param pts
#' @param upper_left
#' @param horizon_radius
#' @param zenith_col_row
#'
#' @returns
#' @export
#'
#' @examples
crop_points <- function(pts,
                         upper_left = NULL,
                         horizon_radius = NULL,
                         zenith_col_row = NULL) {

  .check_vector(upper_left, "numeric", 2, sign = "positive", allow_null = TRUE)

  if (is.null(upper_left)) {
    .check_vector(horizon_radius, "numeric", 1, sign = "positive")
    .check_vector(zenith_col_row, "numeric", 2, sign = "positive")

    pts$col <- pts$col - zenith_col_row[1] + horizon_radius
    pts$row <- pts$row - zenith_col_row[2] + horizon_radius
  } else {
    pts$col <- pts$col - upper_left[1]
    pts$row <- pts$row - upper_left[2]
  }
  pts
}
