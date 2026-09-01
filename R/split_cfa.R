#' Split color filter array
#'
#' @description
#' Split mosaiced sensor data into four layers.
#'
#' @details
#' The input `sensor_data` is expected to contain a single mosaiced sensor
#' layer. The `cfa_pattern` argument defines the channel name associated with
#' each position in the 2 x 2 CFA block:
#' \deqn{
#' \begin{matrix}
#' 0 & 1 \\
#' 2 & 3
#' \end{matrix}
#' }
#'
#' If `black_level` has length one, the same value is subtracted from all CFA
#' channels. If it has length four, values are matched to CFA labels `0`, `1`,
#' `2`, and `3`.
#'
#' @param sensor_data numeric [terra::SpatRaster-class] with one layer.
#' @param cfa_pattern character matrix with two rows and two columns. CFA
#'   channel names associated with labels `0`, `1`, `2`, and `3`. The physical
#'   CFA may contain repeated colors, but each cell must use a unique name
#'   (e.g. `"Green1"`, `"Green2"`).
#' @param cfa_index optional integer [terra::SpatRaster-class] with one layer and the
#'   same dimensions as `sensor_data`. Values must be `0`, `1`, `2`, and `3`.
#'   If `NULL`, labels are generated with [cfa_index_image()].
#' @param black_level optional numeric vector of length one or four. Black level to
#'   subtract before splitting CFA channels. Values must be nonnegative.
#'
#' @return Numeric [terra::SpatRaster-class] with one layer per CFA channel.
#'   Layer names are taken from `cfa_pattern`. The output has half the number
#'   of rows and columns of `sensor_data`.
#'
#' @export
split_cfa <- function(sensor_data, cfa_pattern, cfa_index = NULL, black_level = NULL) {
  .assert_single_layer(sensor_data, allow_na = FALSE)
  .check_cfa_pattern(cfa_pattern)
  if (is.null(cfa_index)) {
    cfa_index <- cfa_index_image(sensor_data)
  } else {
    .assert_single_layer(cfa_index)
    .assert_raster_range(cfa_index, 0, 3)
    .assert_same_geom(sensor_data, cfa_index)
  }
  .check_vector(black_level, "numeric", sign = "nonnegative", allow_null = TRUE)

  raw_pattern <- matrix(c(0,1,
                          2,3), byrow = TRUE, ncol = 2)

  if(is.numeric(black_level)) {
    if (length(black_level) == 1) {
      black_level <- rep(black_level, 4)
    } else if (length(black_level) != 4){
      stop("`black_level` length must be equal to 1 or 4.")
    }
  } else {
    black_level <- c(0, 0, 0, 0)
  }

  caim <- lapply(as.vector(raw_pattern), function(label) {
    sensor_data[cfa_index == label] <-
      sensor_data[cfa_index == label] - black_level[label + 1]
    r <- rast(sensor_data) %>% terra::aggregate(., 2)
    r[] <- sensor_data[cfa_index == label]
    set_rcaiman_geometry(r)
  })
  names(caim) <- as.vector(cfa_pattern)
  rast(caim)
}
